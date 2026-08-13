# ==============================================================================
# Tumor-informed ctDNA across the neoadjuvant continuum in early breast cancer
# Corrected, reproducible meta-analysis and publication-figure workflow
#
# Data source: the accompanying extraction workbook (Effect_Estimates sheet).
# Standard exposure direction: unfavourable ctDNA result vs favourable result.
# Primary statistical model: REML random effects with Hartung-Knapp inference.
# Repeated timepoints: handled with a multilevel model plus CR2 robust inference.
# ==============================================================================

# ------------------------------ 0. Configuration ------------------------------

options(stringsAsFactors = FALSE, scipen = 999)
set.seed(20260801)

INSTALL_MISSING <- TRUE
RUN_BAYESIAN <- TRUE
ASSUMED_WITHIN_STUDY_RHO <- 0.60

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/"))
} else {
  normalizePath(getwd(), winslash = "/")
}

documents_dir <- file.path(Sys.getenv("USERPROFILE"), "Documents")
workbook_candidates <- c(
  file.path(documents_dir, "tumor_informed_ctDNA_meta_extraction_14_studies.xlsx"),
  file.path(script_dir, "tumor_informed_ctDNA_meta_extraction_14_studies.xlsx")
)
workbook_file <- workbook_candidates[file.exists(workbook_candidates)][1]
if (length(workbook_file) == 0 || is.na(workbook_file)) {
  workbook_file <- workbook_candidates[1]
}
analysis_root <- dirname(workbook_file)
out_dir <- file.path(analysis_root, "ctDNA_meta_outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Use a project-local R library. This avoids Windows permission errors when R is
# installed under Program Files and also keeps the analysis environment portable.
r_major_minor <- paste(
  R.version$major,
  strsplit(R.version$minor, "\\.")[[1]][1],
  sep = "."
)
project_library <- file.path(analysis_root, "R_library", r_major_minor)
dir.create(project_library, recursive = TRUE, showWarnings = FALSE)
.libPaths(unique(c(project_library, .libPaths())))

if (identical(getOption("repos")["CRAN"], "@CRAN@") ||
    is.na(getOption("repos")["CRAN"])) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

required_packages <- c(
  "readxl", "dplyr", "tidyr", "purrr", "tibble", "ggplot2", "ggrepel",
  "metafor", "clubSandwich", "svglite", "ragg", "scales"
)
optional_packages <- c(
  "DiagrammeR", "DiagrammeRsvg", "rsvg", "brms", "posterior"
)

install_or_stop <- function(pkgs, required = TRUE, auto_install = required) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) && INSTALL_MISSING && auto_install) {
    message("Installing missing packages into: ", project_library)
    message("Packages: ", paste(missing, collapse = ", "))
    
    # A failed parallel installation can leave 00LOCK folders and cause every
    # dependent package to fail. Only project-local stale locks are removed.
    stale_locks <- list.files(
      project_library, pattern = "^00LOCK", full.names = TRUE,
      recursive = FALSE, all.files = TRUE
    )
    if (length(stale_locks)) {
      message("Removing stale project-library locks: ", paste(basename(stale_locks), collapse = ", "))
      unlink(stale_locks, recursive = TRUE, force = TRUE)
    }
    
    install_log <- file.path(script_dir, "package_install.log")
    retry_order <- unique(c(
      # Explicit dependency-safe order for the packages that failed most often
      # in R 4.3/Windows installations.
      "withr", "scales", "mathjaxr", "pbapply", "digest", "sandwich", "metadat",
      "ggplot2", "metafor", "clubSandwich",
      missing
    ))
    
    install_sequentially <- function() {
      log_connection <- file(install_log, open = "at", encoding = "UTF-8")
      sink(log_connection, type = "output", split = TRUE)
      sink(log_connection, type = "message")
      on.exit({
        sink(type = "message")
        sink(type = "output")
        close(log_connection)
      }, add = TRUE)
      
      cat("\n\n===== Package installation attempt: ", format(Sys.time()), " =====\n", sep = "")
      cat("R version: ", R.version.string, "\n", sep = "")
      cat("Library: ", project_library, "\n", sep = "")
      cat("CRAN: ", getOption("repos")["CRAN"], "\n", sep = "")
      
      for (pkg in retry_order) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
          cat("\n----- Installing ", pkg, " -----\n", sep = "")
          try(
            install.packages(
              pkg,
              lib = project_library,
              dependencies = NA,
              Ncpus = 1L
            ),
            silent = FALSE
          )
        }
      }
    }
    
    install_sequentially()
  }
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) && required) {
    stop(
      "Missing required packages: ", paste(missing, collapse = ", "),
      "\nThe sequential installation retry did not complete.",
      "\n1. Restart R/RStudio and run the script again.",
      "\n2. If it still fails, inspect: ", file.path(script_dir, "package_install.log"),
      "\n3. R 4.3.1 is old; updating R is recommended if CRAN reports version incompatibility."
    )
  }
  invisible(missing)
}

install_or_stop(required_packages, required = TRUE, auto_install = TRUE)
missing_optional <- install_or_stop(optional_packages, required = FALSE, auto_install = FALSE)

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(metafor)
  library(clubSandwich)
  library(scales)
})

if (!file.exists(workbook_file)) {
  stop("Extraction workbook not found: ", workbook_file)
}

# ------------------------------ 1. Figure contract ----------------------------

# Core conclusion:
# Unfavourable tumor-informed ctDNA status is associated with worse time-to-event
# outcomes, with prognostic strength depending on sampling timepoint.
#
# Evidence chain:
# 1) timepoint-specific recurrence-like outcomes;
# 2) adjusted-estimate sensitivity analysis;
# 3) OS, subtype and assay-platform exploratory analyses;
# 4) influence, small-study-effect and Bayesian robustness analyses.
#
# Archetype: quantitative grid, with the timepoint forest as the hero panel.
# Export: editable SVG/PDF plus 600-dpi TIFF and 300-dpi PNG preview.

palette_ctdna <- c(
  baseline = "#3B6FB6",
  during = "#56A6A6",
  preop = "#8A5AA5",
  landmark = "#D98E3D",
  surveillance = "#C95858",
  neutral = "#555555",
  pooled = "#111111"
)

theme_nature <- function(base_size = 7, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.subtitle = element_text(size = base_size, colour = "#444444"),
      strip.text = element_text(face = "bold", size = base_size),
      strip.background = element_rect(fill = "#F3F3F3", colour = NA),
      legend.key.height = grid::unit(3.5, "mm"),
      panel.grid = element_blank(),
      plot.margin = margin(5, 6, 5, 6)
    )
}

theme_set(theme_nature())

# Write PDFs through a temporary file so that an open/locked previous PDF does
# not abort the analysis. Cairo is preferred; the standard PDF device is used
# automatically when Cairo is unavailable on the current Windows/R setup.
safe_pdf_export <- function(final_file, width, height, draw_fun) {
  dir.create(dirname(final_file), recursive = TRUE, showWarnings = FALSE)
  tmp_file <- tempfile(
    pattern = paste0(tools::file_path_sans_ext(basename(final_file)), "_"),
    tmpdir = dirname(final_file),
    fileext = ".pdf"
  )
  on.exit({
    if (file.exists(tmp_file)) unlink(tmp_file, force = TRUE)
  }, add = TRUE)
  
  device_id <- NA_integer_
  cairo_message <- NULL
  tryCatch({
    grDevices::cairo_pdf(
      tmp_file, width = width, height = height,
      family = "sans", onefile = TRUE
    )
    device_id <- grDevices::dev.cur()
  }, error = function(e) {
    cairo_message <<- conditionMessage(e)
  })
  
  if (is.na(device_id)) {
    grDevices::pdf(
      tmp_file, width = width, height = height,
      family = "Helvetica", onefile = TRUE,
      useDingbats = FALSE, compress = TRUE
    )
    device_id <- grDevices::dev.cur()
    warning(
      "Cairo PDF was unavailable; used the standard PDF device. Cairo message: ",
      cairo_message,
      call. = FALSE
    )
  }
  
  draw_error <- tryCatch({
    draw_fun()
    NULL
  }, error = identity)
  
  open_devices <- grDevices::dev.list()
  if (!is.null(open_devices) && device_id %in% open_devices) {
    grDevices::dev.off(device_id)
  }
  if (!is.null(draw_error)) {
    stop(conditionMessage(draw_error), call. = FALSE)
  }
  
  copied <- suppressWarnings(file.copy(tmp_file, final_file, overwrite = TRUE))
  actual_file <- final_file
  if (!isTRUE(copied)) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    actual_file <- sub(
      "\\.pdf$", paste0("_", timestamp, ".pdf"), final_file,
      ignore.case = TRUE
    )
    copied <- suppressWarnings(file.copy(tmp_file, actual_file, overwrite = FALSE))
    if (!isTRUE(copied)) {
      stop(
        "PDF was created but could not be saved in: ", dirname(final_file),
        ". Close open PDF viewers and check folder permissions.",
        call. = FALSE
      )
    }
    warning(
      "The existing PDF appears to be open or locked. Saved instead to: ",
      actual_file,
      call. = FALSE
    )
  }
  
  invisible(actual_file)
}

save_pub <- function(plot, stem, width_mm = 183, height_mm = 120, dpi = 600) {
  stem <- file.path(out_dir, stem)
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  
  svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in)
  print(plot)
  dev.off()
  
  safe_pdf_export(
    paste0(stem, ".pdf"), width_in, height_in,
    function() print(plot)
  )
  
  ragg::agg_tiff(
    paste0(stem, ".tiff"), width = width_in, height = height_in,
    units = "in", res = dpi, compression = "lzw"
  )
  print(plot)
  dev.off()
  
  ragg::agg_png(
    paste0(stem, ".png"), width = width_in, height = height_in,
    units = "in", res = 300
  )
  print(plot)
  dev.off()
}

save_base_plot <- function(plot_fun, stem, width_mm = 120, height_mm = 105,
                           dpi = 600) {
  stem <- file.path(out_dir, stem)
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  
  svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in)
  plot_fun()
  dev.off()
  
  safe_pdf_export(
    paste0(stem, ".pdf"), width_in, height_in,
    plot_fun
  )
  
  ragg::agg_tiff(
    paste0(stem, ".tiff"), width = width_in, height = height_in,
    units = "in", res = dpi, compression = "lzw"
  )
  plot_fun()
  dev.off()
  
  ragg::agg_png(
    paste0(stem, ".png"), width = width_in, height = height_in,
    units = "in", res = 300
  )
  plot_fun()
  dev.off()
}

# Arrange two ggplot objects without patchwork/gridExtra. This uses only R's
# built-in grid package and avoids version conflicts with older ggplot2 builds.
save_pub_pair <- function(plot_left, plot_right, stem, width_mm = 183,
                          height_mm = 95, dpi = 600) {
  stem <- file.path(out_dir, stem)
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  
  draw_pair <- function() {
    grid::grid.newpage()
    layout <- grid::grid.layout(nrow = 1, ncol = 2, widths = grid::unit(c(1, 1), "null"))
    grid::pushViewport(grid::viewport(layout = layout))
    print(
      plot_left,
      vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1),
      newpage = FALSE
    )
    print(
      plot_right,
      vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2),
      newpage = FALSE
    )
    grid::popViewport()
  }
  
  svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in)
  draw_pair()
  dev.off()
  
  safe_pdf_export(
    paste0(stem, ".pdf"), width_in, height_in,
    draw_pair
  )
  
  ragg::agg_tiff(
    paste0(stem, ".tiff"), width = width_in, height = height_in,
    units = "in", res = dpi, compression = "lzw"
  )
  draw_pair()
  dev.off()
  
  ragg::agg_png(
    paste0(stem, ".png"), width = width_in, height = height_in,
    units = "in", res = 300
  )
  draw_pair()
  dev.off()
}

# ------------------------------ 2. Read and standardize data ------------------

raw_effects <- readxl::read_excel(workbook_file, sheet = "Effect_Estimates")
supplementary_results <- readxl::read_excel(workbook_file, sheet = "Supplementary_Results")

effects <- raw_effects %>%
  transmute(
    effect_id = Effect_ID,
    study_id = Study_ID,
    study = Citation,
    cohort = Cohort,
    source_timepoint = Timepoint,
    analysis_n = suppressWarnings(as.numeric(`Analysis N`)),
    subgroup_source = Subgroup,
    endpoint = Endpoint,
    analysis = Analysis,
    reported_contrast = `Reported contrast`,
    hr_reported = as.numeric(HR_reported),
    ci_low_reported = as.numeric(CI_low_reported),
    ci_high_reported = as.numeric(CI_high_reported),
    invert = trimws(as.character(`Invert?`)),
    dependency_cluster = `Dependency cluster`,
    source_locator = `Source locator`,
    notes = Notes
  ) %>%
  mutate(
    hr = if_else(tolower(invert) == "yes", 1 / hr_reported, hr_reported),
    ci_low = if_else(tolower(invert) == "yes", 1 / ci_high_reported, ci_low_reported),
    ci_high = if_else(tolower(invert) == "yes", 1 / ci_low_reported, ci_high_reported),
    yi = log(hr),
    sei = (log(ci_high) - log(ci_low)) / (2 * 1.96),
    vi = sei^2,
    analysis_class = case_when(
      grepl("multivariable|adjusted", analysis, ignore.case = TRUE) ~ "Adjusted",
      TRUE ~ "Unadjusted"
    )
  )

study_dictionary <- tribble(
  ~study_id, ~year, ~assay_family, ~cohort_subtype,
  "S01", 2020, "dPCR-based", "TNBC",
  "S02", 2025, "Sequencing-based", "Mixed",
  "S03", 2025, "Sequencing-based", "Mixed",
  "S04", 2025, "Sequencing-based", "TNBC",
  "S05", 2025, "Sequencing-based", "TNBC",
  "S06", 2026, "Sequencing-based", "Mixed",
  "S07", 2026, "Sequencing-based", "TNBC",
  "S08", 2026, "dPCR-based", "Mixed",
  "S09", 2026, "dPCR-based", "Mixed/non-pCR",
  "S10", 2026, "dPCR-based", "TNBC/non-pCR",
  "S11", 2026, "Sequencing-based", "TNBC + HER2+",
  "S12", 2026, "Sequencing-based", "Mixed",
  "S13", 2025, "Sequencing-based", "Mixed",
  "S14", 2022, "Sequencing-based", "Mixed"
)

effects <- effects %>% left_join(study_dictionary, by = "study_id")

timepoint_map <- tribble(
  ~effect_id, ~clinical_timepoint,
  "E006", "Baseline", "E012", "Baseline", "E061", "Baseline",
  "E062", "Baseline", "E063", "Baseline", "E070", "Baseline",
  "E074", "Baseline",
  "E009", "During NAT", "E023", "During NAT", "E024", "During NAT",
  "E038", "During NAT", "E039", "During NAT", "E064", "During NAT",
  "E065", "During NAT", "E066", "During NAT", "E079", "During NAT",
  "E080", "During NAT", "E083", "During NAT", "E084", "During NAT",
  "E003", "Post-NAT/preoperative", "E004", "Post-NAT/preoperative",
  "E013", "Post-NAT/preoperative", "E014", "Post-NAT/preoperative",
  "E028", "Post-NAT/preoperative", "E036", "Post-NAT/preoperative",
  "E037", "Post-NAT/preoperative", "E050", "Post-NAT/preoperative",
  "E051", "Post-NAT/preoperative", "E052", "Post-NAT/preoperative",
  "E055", "Post-NAT/preoperative", "E056", "Post-NAT/preoperative",
  "E057", "Post-NAT/preoperative", "E058", "Post-NAT/preoperative",
  "E067", "Post-NAT/preoperative", "E068", "Post-NAT/preoperative",
  "E069", "Post-NAT/preoperative", "E071", "Post-NAT/preoperative",
  "E075", "Post-NAT/preoperative", "E076", "Post-NAT/preoperative",
  "E015", "Post-surgery landmark", "E016", "Post-surgery landmark",
  "E040", "Post-surgery landmark", "E041", "Post-surgery landmark",
  "E047", "Post-surgery landmark", "E048", "Post-surgery landmark",
  "E049", "Post-surgery landmark", "E053", "Post-surgery landmark",
  "E054", "Post-surgery landmark", "E059", "Post-surgery landmark",
  "E060", "Post-surgery landmark", "E081", "Post-surgery landmark",
  "E082", "Post-surgery landmark",
  "E005", "Longitudinal surveillance", "E011", "Longitudinal surveillance",
  "E027", "Longitudinal surveillance", "E042", "Longitudinal surveillance",
  "E043", "Longitudinal surveillance", "E044", "Longitudinal surveillance",
  "E077", "Longitudinal surveillance", "E078", "Longitudinal surveillance",
  "E085", "Longitudinal surveillance", "E086", "Longitudinal surveillance"
)

timepoint_levels <- c(
  "Baseline", "During NAT", "Post-NAT/preoperative",
  "Post-surgery landmark", "Longitudinal surveillance"
)

effects <- effects %>%
  left_join(timepoint_map, by = "effect_id") %>%
  mutate(clinical_timepoint = factor(clinical_timepoint, levels = timepoint_levels))

# Analysis sets are prespecified by effect ID. They deliberately exclude:
# - Li 2025 I-SPY2 external validation (overlaps Magbanua 2025 I-SPY2);
# - continuous/burden effects when binary detection is the target exposure;
# - duplicate endpoints from the same cohort and timepoint;
# - estimates without a calculable 95% CI.

primary_timepoint_ids <- c(
  # Baseline
  "E006", "E012", "E062", "E070", "E074",
  # During NAT
  "E009", "E023", "E038", "E065",
  # Post-NAT/preoperative
  "E003", "E014", "E028", "E036", "E052", "E055", "E068", "E071", "E075",
  # Post-surgery landmark
  "E016", "E040", "E047", "E053", "E081",
  # Longitudinal surveillance
  "E005", "E011", "E027", "E042", "E044", "E077", "E086"
)

adjusted_all_ids <- c("E047", "E052", "E055", "E070", "E071", "E076", "E078", "E086")

# One independent estimate per cohort for conventional diagnostics/publication bias.
global_unadjusted_ids <- c(
  "E003", "E005", "E011", "E014", "E081", "E027",
  "E028", "E036", "E044", "E050", "E068", "E075"
)
global_adjusted_ids <- c("E047", "E052", "E055", "E071", "E078", "E086")

# One estimate per each of the 14 independent study cohorts; mixed adjusted/unadjusted.
# This set is exploratory and used only for assay comparison and descriptive diagnostics.
global_all_studies_ids <- c(
  "E003", "E005", "E011", "E014", "E081", "E027", "E028",
  "E036", "E047", "E052", "E055", "E068", "E071", "E075"
)

os_primary_ids <- c(
  "E063",                         # baseline
  "E024", "E039", "E066",      # during NAT
  "E004", "E037", "E051", "E069", # preoperative
  "E041", "E049", "E054", "E082", # post-surgery landmark
  "E043"                          # surveillance
)

subtype_ids <- c(
  "E083",                         # HR+/HER2- during NAT
  "E084", "E023",                # TNBC during NAT
  "E003", "E014", "E028", "E052", "E057", # TNBC preoperative
  "E053", "E059", "E081",      # TNBC postoperative
  "E058", "E060"                 # HER2+ pre/postoperative
)

subtype_labels <- tribble(
  ~effect_id, ~subtype_label,
  "E083", "HR+/HER2-", "E084", "TNBC", "E023", "TNBC",
  "E003", "TNBC", "E014", "TNBC", "E028", "TNBC",
  "E052", "TNBC", "E057", "TNBC", "E053", "TNBC",
  "E059", "TNBC", "E081", "TNBC", "E058", "HER2+", "E060", "HER2+"
)

select_effects <- function(ids) {
  missing_ids <- setdiff(ids, effects$effect_id)
  if (length(missing_ids)) stop("Effect IDs missing from workbook: ", paste(missing_ids, collapse = ", "))
  effects %>% filter(effect_id %in% ids) %>% match_order(ids)
}

match_order <- function(data, ids) {
  data %>% mutate(.order = match(effect_id, ids)) %>% arrange(.order) %>% select(-.order)
}

primary_timepoint <- select_effects(primary_timepoint_ids)
adjusted_all <- select_effects(adjusted_all_ids)
global_unadjusted <- select_effects(global_unadjusted_ids)
global_adjusted <- select_effects(global_adjusted_ids)
global_all_studies <- select_effects(global_all_studies_ids)
os_primary <- select_effects(os_primary_ids)
subtype_data <- select_effects(subtype_ids) %>% left_join(subtype_labels, by = "effect_id")
magbanua_clearance <- select_effects(c("E087", "E088", "E072", "E073")) %>%
  mutate(
    panel = "Full-cohort adjusted clearance trajectory",
    study = reported_contrast
  )
magbanua_rcb <- select_effects(c("E089", "E090", "E091", "E092")) %>%
  mutate(
    panel = if_else(grepl("Pretreatment", source_timepoint), "Pretreatment T0", "Post-NAT/preoperative T3"),
    study = paste0(subgroup_source, ": ctDNA-positive vs negative")
  )

# ------------------------------ 3. Data integrity checks -----------------------

qc_effects <- function(d, name) {
  if (anyDuplicated(d$effect_id)) stop(name, ": duplicated Effect_ID")
  if (any(!is.finite(d$hr) | !is.finite(d$ci_low) | !is.finite(d$ci_high))) {
    stop(name, ": non-finite HR/CI")
  }
  if (any(d$hr <= 0 | d$ci_low <= 0 | d$ci_high <= 0)) stop(name, ": non-positive HR/CI")
  if (any(d$ci_low >= d$ci_high)) stop(name, ": CI lower bound is not below upper bound")
  if (any(d$hr < d$ci_low | d$hr > d$ci_high)) stop(name, ": HR lies outside its CI")
  invisible(TRUE)
}

walk2(
  list(primary_timepoint, adjusted_all, global_unadjusted, global_adjusted,
       global_all_studies, os_primary, subtype_data, magbanua_clearance, magbanua_rcb),
  c("primary_timepoint", "adjusted_all", "global_unadjusted", "global_adjusted",
    "global_all_studies", "os_primary", "subtype_data", "magbanua_clearance", "magbanua_rcb"),
  qc_effects
)

if (anyDuplicated(global_unadjusted$study_id)) stop("global_unadjusted is not independent by cohort")
if (anyDuplicated(global_adjusted$study_id)) stop("global_adjusted is not independent by cohort")
if (anyDuplicated(global_all_studies$study_id)) stop("global_all_studies is not independent by cohort")

dup_primary <- primary_timepoint %>% count(study_id, clinical_timepoint) %>% filter(n > 1)
if (nrow(dup_primary)) stop("Primary timepoint set has duplicate cohort-timepoint rows")

write.csv(effects, file.path(out_dir, "Source_Data_all_standardized_effects.csv"), row.names = FALSE)
write.csv(primary_timepoint, file.path(out_dir, "Source_Data_primary_timepoint.csv"), row.names = FALSE)
write.csv(global_unadjusted, file.path(out_dir, "Source_Data_global_unadjusted.csv"), row.names = FALSE)
write.csv(global_adjusted, file.path(out_dir, "Source_Data_global_adjusted.csv"), row.names = FALSE)
write.csv(os_primary, file.path(out_dir, "Source_Data_OS.csv"), row.names = FALSE)
write.csv(subtype_data, file.path(out_dir, "Source_Data_subtype.csv"), row.names = FALSE)
write.csv(magbanua_clearance,
          file.path(out_dir, "Source_Data_Magbanua_clearance_trajectory.csv"), row.names = FALSE)
write.csv(magbanua_rcb,
          file.path(out_dir, "Source_Data_Magbanua_RCB_strata.csv"), row.names = FALSE)
write.csv(supplementary_results,
          file.path(out_dir, "Source_Data_nonpooled_supplementary_results.csv"), row.names = FALSE)

# ------------------------------ 4. Statistical helpers ------------------------

fit_reml <- function(d, slab = d$study) {
  if (nrow(d) < 2) return(NULL)
  metafor::rma.uni(
    yi = d$yi, vi = d$vi, slab = slab,
    method = "REML", test = "knha"
  )
}

model_row <- function(fit, model_name, group = "Overall") {
  if (is.null(fit)) return(tibble())
  tibble(
    model = model_name,
    group = group,
    k = fit$k,
    log_hr = as.numeric(fit$b),
    hr = exp(as.numeric(fit$b)),
    ci_low = exp(fit$ci.lb),
    ci_high = exp(fit$ci.ub),
    p_value = fit$pval,
    tau2 = fit$tau2,
    i2 = fit$I2,
    q_p_value = fit$QEp,
    method = "REML + Hartung-Knapp"
  )
}

fit_by_group <- function(d, group_col, model_name) {
  group_sym <- rlang::ensym(group_col)
  d %>%
    filter(!is.na(!!group_sym)) %>%
    group_split(!!group_sym, .keep = TRUE) %>%
    map_dfr(function(g) {
      label <- as.character(g %>% pull(!!group_sym) %>% first())
      model_row(fit_reml(g), model_name, label)
    })
}

make_sampling_V <- function(d, rho = 0.60) {
  V <- diag(d$vi)
  for (s in unique(d$study_id)) {
    idx <- which(d$study_id == s)
    if (length(idx) > 1) {
      V[idx, idx] <- rho * sqrt(outer(d$vi[idx], d$vi[idx]))
      diag(V)[idx] <- d$vi[idx]
    }
  }
  V
}

fit_continuum_rve <- function(d, rho = 0.60, moderator = TRUE) {
  d <- d %>% arrange(study_id, clinical_timepoint, effect_id)
  V <- make_sampling_V(d, rho)
  if (moderator) {
    # The model matrix already contains one coefficient per timepoint, so the
    # default intercept must be disabled to avoid redundant predictors.
    mods <- model.matrix(~ 0 + clinical_timepoint, data = d)
    fit <- metafor::rma.mv(
      yi = d$yi, V = V, mods = mods, intercept = FALSE,
      random = ~ 1 | study_id/effect_id,
      method = "REML", data = d
    )
  } else {
    # rma.mv does not accept an explicitly supplied mods=NULL object. Omit the
    # argument entirely for the intercept-only global model.
    fit <- metafor::rma.mv(
      yi = d$yi, V = V,
      random = ~ 1 | study_id/effect_id,
      method = "REML", data = d
    )
  }
  robust <- clubSandwich::coef_test(fit, vcov = "CR2", cluster = d$study_id)
  list(data = d, V = V, fit = fit, robust = robust, rho = rho)
}

continuum_rve <- fit_continuum_rve(primary_timepoint, ASSUMED_WITHIN_STUDY_RHO, TRUE)
continuum_global_rve <- fit_continuum_rve(primary_timepoint, ASSUMED_WITHIN_STUDY_RHO, FALSE)

rho_sensitivity <- map_dfr(c(0, 0.30, 0.60, 0.90), function(rho) {
  obj <- fit_continuum_rve(primary_timepoint, rho, FALSE)
  ct <- as.data.frame(obj$robust)
  critical_t <- qt(0.975, df = ct$df_Satt[1])
  tibble(
    rho = rho,
    estimate_log_hr = ct$beta[1],
    robust_se = ct$SE[1],
    df = ct$df_Satt[1],
    hr = exp(ct$beta[1]),
    ci_low = exp(ct$beta[1] - critical_t * ct$SE[1]),
    ci_high = exp(ct$beta[1] + critical_t * ct$SE[1]),
    p_value = ct$p_Satt[1]
  )
})

write.csv(as.data.frame(continuum_rve$robust),
          file.path(out_dir, "Model_RVE_timepoint_CR2.csv"), row.names = FALSE)
write.csv(as.data.frame(continuum_global_rve$robust),
          file.path(out_dir, "Model_RVE_global_CR2.csv"), row.names = FALSE)
write.csv(rho_sensitivity,
          file.path(out_dir, "Sensitivity_within_study_rho.csv"), row.names = FALSE)

fit_uni_global <- fit_reml(global_unadjusted)
fit_adj_global <- fit_reml(global_adjusted)

model_results <- bind_rows(
  fit_by_group(primary_timepoint, clinical_timepoint, "Primary recurrence-like by timepoint"),
  fit_by_group(os_primary, clinical_timepoint, "OS by timepoint"),
  fit_by_group(subtype_data, subtype_label, "Molecular subtype exploratory"),
  fit_by_group(global_all_studies, assay_family, "Assay platform exploratory"),
  model_row(fit_uni_global, "Independent unadjusted diagnostic set"),
  model_row(fit_adj_global, "Independent adjusted diagnostic set")
)
write.csv(model_results, file.path(out_dir, "Meta_model_results.csv"), row.names = FALSE)

# ------------------------------ 5. Forest-plot helpers ------------------------

format_p_value <- function(p) {
  if (length(p) == 0 || is.na(p) || !is.finite(p)) return("P=NA")
  if (p < 0.001) return("P<0.001")
  sprintf("P=%.3f", p)
}

format_hr_ci_p <- function(hr, ci_low, ci_high, p) {
  sprintf("%.2f (%.2f-%s); %s", hr, ci_low,
          ifelse(ci_high >= 1000, format(ci_high, digits = 3, scientific = TRUE),
                 sprintf("%.2f", ci_high)),
          format_p_value(p))
}

fit_annotation <- function(fit) {
  sprintf("Pooled HR %.2f (95%% CI %.2f-%.2f); %s",
          exp(as.numeric(fit$b)), exp(fit$ci.lb), exp(fit$ci.ub),
          format_p_value(fit$pval))
}

forest_plot <- function(d, group_col, title, subtitle = NULL, pool = TRUE,
                        colour_map = NULL) {
  group_sym <- rlang::ensym(group_col)
  d <- d %>%
    filter(!is.na(!!group_sym)) %>%
    mutate(
      group_plot = as.character(!!group_sym),
      group_plot = factor(group_plot, levels = unique(group_plot))
    )
  
  plot_parts <- d %>%
    group_split(group_plot, .keep = TRUE) %>%
    map(function(g) {
      fit <- fit_reml(g)
      if (!is.null(fit)) {
        w <- as.numeric(weights(fit))
        g$weight_plot <- 1.8 + 3.2 * sqrt(w / max(w))
      } else {
        g$weight_plot <- 3
      }
      g <- g %>%
        mutate(
          is_pooled = FALSE,
          display_label = paste0(study, " (", endpoint, "; ", analysis_class, ")"),
          p_value = 2 * pnorm(-abs(yi / sei)),
          stat_label = purrr::pmap_chr(
            list(hr, ci_low, ci_high, p_value), format_hr_ci_p
          )
        )
      if (pool && !is.null(fit)) {
        pooled <- tibble(
          effect_id = paste0("POOL_", unique(g$group_plot)),
          study_id = NA_character_, study = NA_character_, endpoint = NA_character_,
          analysis_class = NA_character_, group_plot = unique(g$group_plot),
          hr = exp(as.numeric(fit$b)), ci_low = exp(fit$ci.lb), ci_high = exp(fit$ci.ub),
          weight_plot = 5.2, is_pooled = TRUE,
          p_value = as.numeric(fit$pval),
          stat_label = format_hr_ci_p(
            exp(as.numeric(fit$b)), exp(fit$ci.lb), exp(fit$ci.ub), fit$pval
          ),
          display_label = sprintf(
            "Random effects (k=%d; I%s=%.1f%%)", fit$k, "\u00b2", fit$I2
          )
        )
        bind_rows(g, pooled)
      } else {
        g
      }
    }) %>%
    bind_rows() %>%
    group_by(group_plot) %>%
    arrange(is_pooled, study, .by_group = TRUE) %>%
    mutate(row_index = row_number()) %>%
    ungroup() %>%
    mutate(
      row_key = paste(group_plot, row_index, display_label, sep = "___"),
      row_key = factor(row_key, levels = rev(unique(row_key))),
      stat_font = if_else(is_pooled, "bold", "plain")
    )
  
  label_lookup <- setNames(plot_parts$display_label, as.character(plot_parts$row_key))
  subtitle_full <- paste(
    c(subtitle, "Right column: HR (95% CI); P. Individual P values are two-sided Wald values derived from reported HR and 95% CI."),
    collapse = "\n"
  )
  
  p <- ggplot(plot_parts, aes(x = hr, y = row_key)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "#777777", linewidth = 0.4) +
    geom_segment(
      aes(x = ci_low, xend = ci_high, yend = row_key, colour = group_plot),
      linewidth = 0.55
    ) +
    geom_point(
      aes(size = weight_plot, shape = is_pooled, colour = group_plot),
      stroke = 0.35
    ) +
    geom_text(
      aes(x = Inf, label = stat_label, fontface = stat_font),
      hjust = -0.04, size = 2.0, colour = "#222222", show.legend = FALSE
    ) +
    facet_grid(rows = vars(group_plot), scales = "free_y", space = "free_y", switch = "y") +
    scale_x_log10(
      breaks = c(0.25, 1, 5, 25, 100, 500, 2500),
      labels = scales::label_number(accuracy = 0.01)
    ) +
    scale_y_discrete(labels = label_lookup) +
    scale_shape_manual(values = c(`FALSE` = 15, `TRUE` = 18), guide = "none") +
    scale_size_identity() +
    coord_cartesian(clip = "off") +
    labs(x = "Hazard ratio (95% CI)", y = NULL, title = title, subtitle = subtitle_full) +
    theme_nature(base_size = 6.5) +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 1),
      panel.spacing.y = grid::unit(2.5, "mm"),
      legend.position = "none",
      plot.margin = margin(5.5, 68, 5.5, 5.5, unit = "mm")
    )
  
  if (!is.null(colour_map)) {
    p <- p + scale_colour_manual(values = colour_map)
  } else {
    groups <- unique(plot_parts$group_plot)
    auto_cols <- setNames(rep(c("#3B6FB6", "#8A5AA5", "#D98E3D", "#C95858", "#56A6A6"),
                              length.out = length(groups)), groups)
    p <- p + scale_colour_manual(values = auto_cols)
  }
  p
}

timepoint_colours <- c(
  "Baseline" = unname(palette_ctdna["baseline"]),
  "During NAT" = unname(palette_ctdna["during"]),
  "Post-NAT/preoperative" = unname(palette_ctdna["preop"]),
  "Post-surgery landmark" = unname(palette_ctdna["landmark"]),
  "Longitudinal surveillance" = unname(palette_ctdna["surveillance"])
)

# ------------------------------ 6. Main figures -------------------------------

fig1a <- forest_plot(
  primary_timepoint, clinical_timepoint,
  "Recurrence-like outcomes across the treatment continuum",
  "One effect per cohort within each timepoint; REML random effects with Hartung-Knapp inference",
  pool = TRUE, colour_map = timepoint_colours
)
save_pub(fig1a, "Figure_1A_Timepoint_Forest", 230, 220)

fig1b <- forest_plot(
  adjusted_all, clinical_timepoint,
  "Adjusted hazard-ratio estimates",
  "Adjusted models are pooled only within clinically aligned timepoint strata",
  pool = TRUE, colour_map = timepoint_colours
)
save_pub(fig1b, "Figure_1B_Adjusted_Forest", 230, 125)

fig2 <- forest_plot(
  os_primary, clinical_timepoint,
  "Overall survival",
  "No across-timepoint conventional pooled estimate is shown",
  pool = TRUE, colour_map = timepoint_colours
)
save_pub(fig2, "Figure_2_OS_Forest", 230, 145)

subtype_data <- subtype_data %>%
  mutate(subtype_timepoint = paste(subtype_label, clinical_timepoint, sep = " | "))
fig3a <- forest_plot(
  subtype_data, subtype_timepoint,
  "Molecular-subtype analysis",
  "Exploratory: adjusted and unadjusted estimates are identified in row labels",
  pool = TRUE
)
save_pub(fig3a, "Figure_3A_Subtype_Forest", 230, 150)

fig3b <- forest_plot(
  global_all_studies, assay_family,
  "Assay-platform subgroup analysis",
  "Exploratory and potentially confounded by timepoint, cohort and analysis type",
  pool = TRUE,
  colour_map = c("Sequencing-based" = "#3B6FB6", "dPCR-based" = "#D98E3D")
)
save_pub(fig3b, "Figure_3B_Assay_Forest", 230, 125)

# ------------------------------ 7. Funnel and small-study effects -------------

funnel_plot <- function(d, fit, title, egger_p = NA_real_) {
  mu <- as.numeric(fit$b)
  max_se <- max(d$sei) * 1.06
  boundary <- tibble(
    sei = seq(0, max_se, length.out = 250),
    left = exp(mu - 1.96 * sei),
    right = exp(mu + 1.96 * sei)
  )
  # Construct the pseudo-95% region as a closed polygon. This is compatible
  # with older ggplot2 versions that cannot draw a horizontal geom_ribbon.
  funnel_region <- bind_rows(
    boundary %>% transmute(x = left, y = sei),
    boundary %>% arrange(desc(sei)) %>% transmute(x = right, y = sei)
  )
  pooled_text <- fit_annotation(fit)
  subtitle <- if (is.finite(egger_p)) {
    sprintf("%s\nEgger intercept %s; k=%d", pooled_text,
            format_p_value(egger_p), nrow(d))
  } else {
    sprintf("%s\nEgger test not performed because k=%d (<10)",
            pooled_text, nrow(d))
  }
  ggplot(d, aes(x = hr, y = sei)) +
    geom_polygon(
      data = funnel_region,
      aes(x = x, y = y), inherit.aes = FALSE,
      fill = "#DDE8F2", colour = NA, alpha = 0.65
    ) +
    geom_vline(xintercept = exp(mu), colour = "#333333", linewidth = 0.55) +
    geom_line(data = boundary, aes(x = left, y = sei), inherit.aes = FALSE,
              linetype = 2, colour = "#777777") +
    geom_line(data = boundary, aes(x = right, y = sei), inherit.aes = FALSE,
              linetype = 2, colour = "#777777") +
    geom_point(shape = 21, size = 2.5, fill = "#3B6FB6", colour = "black", stroke = 0.3) +
    ggrepel::geom_text_repel(aes(label = study), size = 2.0, max.overlaps = Inf,
                             min.segment.length = 0, seed = 20260801) +
    scale_x_log10() +
    scale_y_reverse(expand = expansion(mult = c(0.03, 0.08))) +
    labs(x = "Hazard ratio (log scale)", y = "Standard error", title = title, subtitle = subtitle) +
    theme_nature(base_size = 7)
}

egger_standard <- function(d) {
  if (nrow(d) < 10) return(NULL)
  egger_df <- d %>% mutate(precision = 1 / sei, snd = yi / sei)
  fit <- lm(snd ~ precision, data = egger_df)
  list(
    data = egger_df,
    fit = fit,
    intercept = unname(coef(fit)[1]),
    p_value = coef(summary(fit))[1, "Pr(>|t|)"]
  )
}

egger_uni <- egger_standard(global_unadjusted)
egger_adj <- egger_standard(global_adjusted)

fig4a <- funnel_plot(
  global_unadjusted, fit_uni_global, "Funnel plot: independent unadjusted estimates",
  if (is.null(egger_uni)) NA_real_ else egger_uni$p_value
)
save_pub(fig4a, "Figure_4A_Funnel_Unadjusted", 120, 105)

fig4b <- funnel_plot(
  global_adjusted, fit_adj_global, "Funnel plot: independent adjusted estimates",
  if (is.null(egger_adj)) NA_real_ else egger_adj$p_value
)
save_pub(fig4b, "Figure_4B_Funnel_Adjusted", 120, 105)

# Proper Egger plot: standardized normal deviate vs precision.
if (!is.null(egger_uni)) {
  egger_line <- tibble(
    precision = seq(min(egger_uni$data$precision), max(egger_uni$data$precision), length.out = 200)
  )
  pred <- predict(egger_uni$fit, newdata = egger_line, interval = "confidence")
  egger_line <- bind_cols(egger_line, as.data.frame(pred))
  
  p_egger <- ggplot(egger_uni$data, aes(x = precision, y = snd)) +
    geom_ribbon(
      data = egger_line, aes(x = precision, ymin = lwr, ymax = upr),
      inherit.aes = FALSE,
      fill = "#E7B6B2", alpha = 0.35
    ) +
    geom_line(data = egger_line, aes(y = fit), colour = "#C95858", linewidth = 0.7) +
    geom_point(shape = 21, size = 2.4, fill = "#3B6FB6", colour = "black", stroke = 0.3) +
    ggrepel::geom_text_repel(aes(label = study), size = 2.0, max.overlaps = Inf,
                             min.segment.length = 0, seed = 20260801) +
    geom_hline(yintercept = 0, linetype = 2, colour = "#777777") +
    labs(
      x = "Precision (1/SE)", y = "Standard normal deviate (log HR/SE)",
      title = "Egger regression for funnel-plot asymmetry",
      subtitle = sprintf("Intercept=%.2f (95%% confidence band); %s; k=%d",
                         egger_uni$intercept, format_p_value(egger_uni$p_value),
                         nrow(egger_uni$data))
    ) +
    theme_nature(base_size = 7)
  save_pub(p_egger, "Supplementary_Figure_9_Egger_Regression", 135, 105)
  save_pub(p_egger, "Egger_Regression_Plot_Final", 135, 105)
  capture.output(summary(egger_uni$fit),
                 file = file.path(out_dir, "Egger_regression_model.txt"))
}

# ------------------------------ 8. Cumulative analyses ------------------------

cumulative_data <- function(d) {
  d <- d %>% arrange(year, study)
  map_dfr(seq_len(nrow(d)), function(i) {
    di <- d[seq_len(i), , drop = FALSE]
    if (i == 1) {
      tibble(
        step = i, year = di$year[i], study = di$study[i], k = 1,
        hr = di$hr[i], ci_low = di$ci_low[i], ci_high = di$ci_high[i],
        p_value = 2 * pnorm(-abs(di$yi[i] / di$sei[i]))
      )
    } else {
      fit <- fit_reml(di)
      tibble(
        step = i, year = di$year[i], study = di$study[i], k = i,
        hr = exp(as.numeric(fit$b)), ci_low = exp(fit$ci.lb), ci_high = exp(fit$ci.ub),
        p_value = as.numeric(fit$pval)
      )
    }
  }) %>%
    mutate(
      label = paste0(year, "  ", study),
      stat_label = purrr::pmap_chr(
        list(hr, ci_low, ci_high, p_value), format_hr_ci_p
      ),
      hr_plot = pmin(pmax(hr, 0.25), 500),
      ci_low_plot = pmax(ci_low, 0.25),
      ci_high_plot = pmin(ci_high, 500)
    )
}

plot_cumulative <- function(d, title) {
  d <- d %>% mutate(label = factor(label, levels = rev(label)))
  ggplot(d, aes(x = hr_plot, y = label)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "#777777", linewidth = 0.4) +
    geom_segment(aes(x = ci_low_plot, xend = ci_high_plot, yend = label),
                 colour = "#3B6FB6", linewidth = 0.55) +
    geom_point(shape = 18, size = 2.7, colour = "#222222") +
    geom_text(aes(x = Inf, label = stat_label), hjust = -0.04,
              size = 2.0, colour = "#222222") +
    scale_x_log10(limits = c(0.25, 500),
                  breaks = c(0.25, 1, 5, 25, 100, 500),
                  labels = scales::label_number(accuracy = 0.01)) +
    coord_cartesian(clip = "off") +
    labs(x = "Cumulative hazard ratio (95% CI)", y = NULL, title = title,
         subtitle = paste0(
           "Studies added chronologically; REML + Hartung-Knapp after the first study. ",
           "Right column: HR (95% CI); P. CIs outside 0.25-500 are truncated only in the plotting panel."
         )) +
    theme_nature(base_size = 7) +
    theme(plot.margin = margin(5.5, 66, 5.5, 5.5, unit = "mm"))
}

cum_uni <- cumulative_data(global_unadjusted)
cum_adj <- cumulative_data(global_adjusted)
write.csv(cum_uni, file.path(out_dir, "Cumulative_unadjusted_results.csv"), row.names = FALSE)
write.csv(cum_adj, file.path(out_dir, "Cumulative_adjusted_results.csv"), row.names = FALSE)

save_pub(plot_cumulative(cum_uni, "Cumulative meta-analysis: unadjusted estimates"),
         "Figure_5A_Cumulative_Unadjusted", 215, 120)
save_pub(plot_cumulative(cum_adj, "Cumulative meta-analysis: adjusted estimates"),
         "Supplementary_Figure_6_Cumulative_Adjusted", 210, 95)

# ------------------------------ 9. Influence diagnostics ----------------------

save_base_plot(
  function() {
    par(mar = c(5.0, 4.2, 2.0, 1.0), family = "sans")
    metafor::baujat(fit_uni_global, symbol = "ids", cex = 0.85,
                    xlab = expression(Delta*Q), ylab = expression(Influence~on~pooled~effect))
    title("Baujat plot: unadjusted independent set",
          sub = fit_annotation(fit_uni_global), cex.main = 0.95, cex.sub = 0.70)
  },
  "Figure_6A_Baujat_Unadjusted", 120, 105
)

save_base_plot(
  function() {
    par(mar = c(5.0, 4.2, 2.0, 1.0), family = "sans")
    metafor::baujat(fit_adj_global, symbol = "ids", cex = 0.85,
                    xlab = expression(Delta*Q), ylab = expression(Influence~on~pooled~effect))
    title("Baujat plot: adjusted independent set",
          sub = fit_annotation(fit_adj_global), cex.main = 0.95, cex.sub = 0.70)
  },
  "Figure_6B_Baujat_Adjusted", 120, 105
)

save_base_plot(
  function() {
    par(mar = c(5.0, 4.2, 2.0, 1.0), family = "sans")
    metafor::radial(fit_uni_global)
    title("Radial plot: unadjusted independent set",
          sub = fit_annotation(fit_uni_global), cex.main = 0.95, cex.sub = 0.70)
  },
  "Figure_7A_Radial_Unadjusted", 120, 105
)

save_base_plot(
  function() {
    par(mar = c(5.0, 4.2, 2.0, 1.0), family = "sans")
    metafor::radial(fit_adj_global)
    title("Radial plot: adjusted independent set",
          sub = fit_annotation(fit_adj_global), cex.main = 0.95, cex.sub = 0.70)
  },
  "Figure_7B_Radial_Adjusted", 120, 105
)

leave_one_out <- map_dfr(seq_len(nrow(global_unadjusted)), function(i) {
  di <- global_unadjusted[-i, , drop = FALSE]
  fit <- fit_reml(di)
  tibble(
    omitted = global_unadjusted$study[i],
    hr = exp(as.numeric(fit$b)),
    ci_low = exp(fit$ci.lb),
    ci_high = exp(fit$ci.ub),
    p_value = as.numeric(fit$pval),
    i2 = fit$I2
  )
})
write.csv(leave_one_out, file.path(out_dir, "Leave_one_out_results.csv"), row.names = FALSE)

p_loo <- leave_one_out %>%
  mutate(
    omitted = factor(omitted, levels = rev(omitted)),
    stat_label = purrr::pmap_chr(
      list(hr, ci_low, ci_high, p_value), format_hr_ci_p
    )
  ) %>%
  ggplot(aes(x = hr, y = omitted)) +
  geom_vline(xintercept = exp(as.numeric(fit_uni_global$b)), linetype = 2,
             colour = "#777777", linewidth = 0.4) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = omitted),
               colour = "#3B6FB6", linewidth = 0.55) +
  geom_point(shape = 18, size = 2.7) +
  geom_text(aes(x = Inf, label = stat_label), hjust = -0.04,
            size = 2.0, colour = "#222222") +
  scale_x_log10(breaks = c(1, 2, 5, 10, 25, 50, 100),
                labels = scales::label_number(accuracy = 0.01)) +
  coord_cartesian(clip = "off") +
  labs(x = "Hazard ratio after omission (95% CI)", y = "Omitted cohort",
       title = "Leave-one-cohort-out sensitivity analysis",
       subtitle = "Right column: pooled HR (95% CI); Hartung-Knapp P after each omission") +
  theme_nature(base_size = 7) +
  theme(plot.margin = margin(5.5, 65, 5.5, 5.5, unit = "mm"))
save_pub(p_loo, "Supplementary_Figure_5_Leave_One_Out", 210, 115)

# ------------------------------ 10. Timepoint-specific supplementary forests --

timepoint_file_tags <- c(
  "Baseline" = "2A_Baseline",
  "During NAT" = "2B_During_NAT",
  "Post-NAT/preoperative" = "2C_Preoperative",
  "Post-surgery landmark" = "2D_Postsurgery_Landmark",
  "Longitudinal surveillance" = "2E_Surveillance"
)

for (tp in timepoint_levels) {
  d_tp <- primary_timepoint %>% filter(clinical_timepoint == tp) %>% mutate(panel = tp)
  if (nrow(d_tp)) {
    p_tp <- forest_plot(
      d_tp, panel, paste0(tp, " ctDNA"),
      "Recurrence-like time-to-event outcomes", pool = TRUE,
      colour_map = setNames(timepoint_colours[tp], tp)
    )
    save_pub(p_tp, paste0("Supplementary_Figure_", timepoint_file_tags[tp]), 220,
             max(70, 45 + 8 * nrow(d_tp)))
  }
}

# Standalone subtype panels
for (st in unique(subtype_data$subtype_label)) {
  d_st <- subtype_data %>% filter(subtype_label == st) %>% mutate(panel = st)
  p_st <- forest_plot(
    d_st, panel, paste0(st, " subtype"),
    "Exploratory subtype-specific effects", pool = nrow(d_st) >= 2
  )
  safe_st <- gsub("[^A-Za-z0-9]+", "_", st)
  save_pub(p_st, paste0("Supplementary_Figure_3_", safe_st), 220,
           max(65, 45 + 8 * nrow(d_st)))
}

# Standalone assay panels
for (assay in unique(global_all_studies$assay_family)) {
  d_assay <- global_all_studies %>% filter(assay_family == assay) %>% mutate(panel = assay)
  p_assay <- forest_plot(
    d_assay, panel, assay,
    "Exploratory one-estimate-per-cohort analysis", pool = nrow(d_assay) >= 2
  )
  safe_assay <- gsub("[^A-Za-z0-9]+", "_", assay)
  save_pub(p_assay, paste0("Supplementary_Figure_4_", safe_assay), 220,
           max(65, 45 + 8 * nrow(d_assay)))
}

# Magbanua 2025: correlated estimates from one cohort are displayed but not pooled.
p_mag_clearance <- forest_plot(
  magbanua_clearance, panel,
  "Magbanua 2025: adjusted ctDNA clearance trajectory",
  "Common reference: persistently ctDNA-negative; correlated estimates are not meta-analysed",
  pool = FALSE,
  colour_map = c("Full-cohort adjusted clearance trajectory" = "#8A5AA5")
)
save_pub(p_mag_clearance, "Supplementary_Figure_12_Magbanua_Clearance_Trajectory", 220, 85)

p_mag_rcb <- forest_plot(
  magbanua_rcb, panel,
  "Magbanua 2025: RCB-II/III risk stratification",
  "Single-cohort stratified estimates; no pooled summary is calculated",
  pool = FALSE,
  colour_map = c(
    "Pretreatment T0" = "#3B6FB6",
    "Post-NAT/preoperative T3" = "#8A5AA5"
  )
)
save_pub(p_mag_rcb, "Supplementary_Figure_13_Magbanua_RCB_Strata", 220, 90)

# ------------------------------ 11. Trim-and-fill -----------------------------

tf_uni <- metafor::trimfill(fit_uni_global)

tf_original <- funnel_plot(
  global_unadjusted, fit_uni_global, "Original funnel plot",
  if (is.null(egger_uni)) NA_real_ else egger_uni$p_value
)

tf_data <- tibble(
  study = ifelse(tf_uni$fill, "Imputed study", "Observed study"),
  yi = as.numeric(tf_uni$yi),
  sei = sqrt(as.numeric(tf_uni$vi)),
  hr = exp(yi),
  imputed = as.logical(tf_uni$fill)
)

tf_adjusted <- ggplot(tf_data, aes(x = hr, y = sei, fill = imputed)) +
  geom_vline(xintercept = exp(as.numeric(tf_uni$b)), colour = "#333333", linewidth = 0.55) +
  geom_point(shape = 21, size = 2.6, colour = "black", stroke = 0.3) +
  scale_fill_manual(values = c(`FALSE` = "#3B6FB6", `TRUE` = "white"),
                    labels = c("Observed", "Imputed")) +
  scale_x_log10() +
  scale_y_reverse() +
  labs(x = "Hazard ratio (log scale)", y = "Standard error",
       title = "Trim-and-fill",
       subtitle = sprintf("Imputed studies=%d; %s", tf_uni$k0, fit_annotation(tf_uni)),
       fill = NULL) +
  theme_nature(base_size = 7) +
  theme(legend.position = "bottom")

save_pub_pair(
  tf_original, tf_adjusted,
  "Supplementary_Figure_10_Trim_and_Fill", 230, 95
)

# ------------------------------ 12. Study weights and HR distribution ---------

weight_data <- global_unadjusted %>%
  mutate(
    random_effect_weight = as.numeric(weights(fit_uni_global)),
    study = factor(study, levels = study[order(random_effect_weight)])
  )
write.csv(weight_data, file.path(out_dir, "Random_effect_weights.csv"), row.names = FALSE)

p_weights <- ggplot(weight_data, aes(x = random_effect_weight, y = study,
                                     fill = clinical_timepoint)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = timepoint_colours, drop = FALSE) +
  labs(x = "Random-effects weight (%)", y = NULL,
       title = "Study weights in the independent unadjusted model",
       subtitle = fit_annotation(fit_uni_global), fill = "Timepoint") +
  theme_nature(base_size = 7) +
  theme(legend.position = "bottom")
save_pub(p_weights, "Supplementary_Figure_8_Study_Weights", 155, 110)

rve_display <- rho_sensitivity %>%
  slice(which.min(abs(rho - ASSUMED_WITHIN_STUDY_RHO)))

distribution_subtitle <- sprintf(
  "CR2 global HR %.2f (95%% CI %.2f-%.2f); %s.\nDistribution remains descriptive because repeated cohort estimates are dependent.",
  rve_display$hr, rve_display$ci_low, rve_display$ci_high,
  format_p_value(rve_display$p_value)
)

p_distribution <- ggplot(
  primary_timepoint,
  aes(x = clinical_timepoint, y = hr, colour = clinical_timepoint)
) +
  geom_boxplot(width = 0.45, outlier.shape = NA, colour = "#444444", fill = NA,
               linewidth = 0.45) +
  geom_jitter(width = 0.12, height = 0, size = 2.2, alpha = 0.85) +
  scale_y_log10() +
  scale_x_discrete(labels = c(
    "Baseline" = "Baseline",
    "During NAT" = "During\nNAT",
    "Post-NAT/preoperative" = "Post-NAT/\npreoperative",
    "Post-surgery landmark" = "Post-surgery\nlandmark",
    "Longitudinal surveillance" = "Longitudinal\nsurveillance"
  )) +
  scale_colour_manual(values = timepoint_colours, guide = "none") +
  labs(x = NULL, y = "Hazard ratio (log scale)",
       title = "Distribution of effect estimates by ctDNA timepoint",
       subtitle = distribution_subtitle) +
  theme_nature(base_size = 7) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, lineheight = 0.9))
save_pub(p_distribution, "Supplementary_Figure_11_Timepoint_Distribution", 183, 105)

# ------------------------------ 13. Bayesian sensitivity analysis -------------

if (RUN_BAYESIAN &&
    requireNamespace("brms", quietly = TRUE) &&
    requireNamespace("posterior", quietly = TRUE)) {
  
  bayes_fit <- brms::brm(
    yi | se(sei, sigma = TRUE) ~ 1,
    data = global_unadjusted,
    family = gaussian(),
    prior = c(
      brms::prior(normal(0, 2), class = Intercept),
      brms::prior(exponential(1), class = sigma)
    ),
    chains = 4, iter = 4000, warmup = 2000,
    seed = 20260801, cores = min(4, parallel::detectCores()),
    control = list(adapt_delta = 0.99, max_treedepth = 12),
    refresh = 200
  )
  
  draws <- posterior::as_draws_df(bayes_fit) %>%
    transmute(hr = exp(b_Intercept), tau = sigma)
  
  bayes_summary <- tibble(
    posterior_median_hr = median(draws$hr),
    ci_low = quantile(draws$hr, 0.025),
    ci_high = quantile(draws$hr, 0.975),
    posterior_median_tau = median(draws$tau),
    tau_ci_low = quantile(draws$tau, 0.025),
    tau_ci_high = quantile(draws$tau, 0.975),
    prob_hr_gt_1 = mean(draws$hr > 1)
  )
  write.csv(bayes_summary, file.path(out_dir, "Bayesian_meta_summary.csv"), row.names = FALSE)
  saveRDS(bayes_fit, file.path(out_dir, "Bayesian_meta_model.rds"))
  
  p_bayes <- ggplot(draws, aes(x = hr)) +
    geom_density(fill = "#3B6FB6", colour = "#254B7A", alpha = 0.55, linewidth = 0.6) +
    geom_vline(xintercept = median(draws$hr), linetype = 2, colour = "#C95858", linewidth = 0.65) +
    scale_x_log10() +
    labs(
      x = "Hazard ratio (log scale)", y = "Posterior density",
      title = "Bayesian random-effects sensitivity analysis",
      subtitle = sprintf("Posterior median HR %.2f (95%% CrI %.2f-%.2f); Pr(HR>1)=%.3f",
                         bayes_summary$posterior_median_hr,
                         bayes_summary$ci_low, bayes_summary$ci_high,
                         bayes_summary$prob_hr_gt_1)
    ) +
    theme_nature(base_size = 7)
  save_pub(p_bayes, "Supplementary_Figure_7_Bayesian_Posterior", 135, 100)
} else {
  message(
    "Bayesian analysis skipped. Install brms and posterior, then set RUN_BAYESIAN <- TRUE."
  )
}

# ------------------------------ 14. PRISMA flow diagram -----------------------

# Do not reuse unverified historical screening counts. Enter the final counts below
# or edit the automatically written template CSV, then rerun the script.
prisma_counts <- tibble(
  identified_databases = NA_integer_,
  identified_other = NA_integer_,
  duplicates_removed = NA_integer_,
  records_screened = NA_integer_,
  records_excluded = NA_integer_,
  reports_sought = NA_integer_,
  reports_not_retrieved = NA_integer_,
  reports_assessed = NA_integer_,
  reports_excluded = NA_integer_,
  studies_included = 14L
)
write.csv(prisma_counts, file.path(out_dir, "PRISMA_counts_template.csv"), row.names = FALSE)

if (all(!is.na(prisma_counts)) &&
    all(c("DiagrammeR", "DiagrammeRsvg", "rsvg") %in%
        rownames(installed.packages()))) {
  pc <- prisma_counts[1, ]
  graph_code <- sprintf(
    "digraph prisma {
       graph [layout=dot, rankdir=TB]
       node [shape=box, style=filled, fillcolor='#EAF2F8', fontname='Arial']
       a [label='Records identified\\nDatabases: n=%d; other: n=%d']
       b [label='Duplicates removed\\nn=%d']
       c [label='Records screened\\nn=%d']
       d [label='Records excluded\\nn=%d']
       e [label='Reports sought\\nn=%d']
       f [label='Reports not retrieved\\nn=%d']
       g [label='Full-text reports assessed\\nn=%d']
       h [label='Full-text reports excluded\\nn=%d']
       i [label='Studies included\\nn=%d', fillcolor='#D5F5E3']
       a -> b -> c; c -> d; c -> e; e -> f; e -> g; g -> h; g -> i
     }",
    pc$identified_databases, pc$identified_other, pc$duplicates_removed,
    pc$records_screened, pc$records_excluded, pc$reports_sought,
    pc$reports_not_retrieved, pc$reports_assessed, pc$reports_excluded,
    pc$studies_included
  )
  gr <- DiagrammeR::grViz(graph_code)
  svg_text <- DiagrammeRsvg::export_svg(gr)
  svg_file <- file.path(out_dir, "Supplementary_Figure_1_PRISMA.svg")
  writeLines(svg_text, svg_file, useBytes = TRUE)
  rsvg::rsvg_pdf(charToRaw(svg_text), file.path(out_dir, "Supplementary_Figure_1_PRISMA.pdf"))
  rsvg::rsvg_png(charToRaw(svg_text), file.path(out_dir, "Supplementary_Figure_1_PRISMA.png"),
                 width = 2400)
} else {
  message("PRISMA diagram not generated: final screening counts and/or optional packages are missing.")
}

# ------------------------------ 15. Audit trail and completion ----------------

capture.output(
  list(
    primary_timepoint_models = fit_by_group(
      primary_timepoint, clinical_timepoint, "Primary recurrence-like by timepoint"
    ),
    continuum_RVE_CR2 = continuum_rve$robust,
    rho_sensitivity = rho_sensitivity,
    independent_unadjusted = summary(fit_uni_global),
    independent_adjusted = summary(fit_adj_global),
    trim_and_fill = summary(tf_uni)
  ),
  file = file.path(out_dir, "Full_model_audit.txt")
)

capture.output(sessionInfo(), file = file.path(out_dir, "SessionInfo.txt"))

cat("\nAnalysis completed.\n")
cat("Input workbook: ", workbook_file, "\n", sep = "")
cat("Output directory: ", out_dir, "\n", sep = "")
cat("Important: across-timepoint conventional pooling was replaced by multilevel CR2 analysis.\n")
cat("Important: Egger testing was restricted to an independent set with k >= 10.\n")
cat("Important: PRISMA output requires the final screening counts in the template.\n")

