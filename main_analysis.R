# Tumor-informed ctDNA meta-analysis: reproducible analysis script
#
# Effect-level data are embedded below; no external data file is required.
#
# Optional environment variables:
#   CTDNA_META_OUTPUT_DIR      output directory (default: ./ctDNA_meta_outputs)
#   CTDNA_META_MAKE_FIGURES   true/false (default: true)
#   CTDNA_META_RUN_BAYESIAN   true/false (default: true)
#   CTDNA_META_INSTALL_MISSING true/false (default: false)
#
# The analysis was validated with R 4.3.1. Package versions are written to
# SessionInfo.txt at the end of each run.

# ------------------------------ 0. Configuration ------------------------------
options(stringsAsFactors = FALSE, scipen = 999)
set.seed(20260801)
INSTALL_MISSING <- identical(
  tolower(Sys.getenv("CTDNA_META_INSTALL_MISSING", unset = "false")), "true"
)
RUN_BAYESIAN <- identical(
  tolower(Sys.getenv("CTDNA_META_RUN_BAYESIAN", unset = "true")), "true"
)
MAKE_FIGURES <- identical(
  tolower(Sys.getenv("CTDNA_META_MAKE_FIGURES", unset = "true")), "true"
)
ASSUMED_WITHIN_COHORT_RHO <- 0.60
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/"))
} else {
  normalizePath(getwd(), winslash = "/")
}
output_from_env <- Sys.getenv("CTDNA_META_OUTPUT_DIR", unset = "")
out_dir <- if (nzchar(output_from_env)) output_from_env else file.path(script_dir, "ctDNA_meta_outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# Use a project-local R library. This avoids Windows permission errors when R is
# installed under Program Files and also keeps the analysis environment portable.
r_major_minor <- paste(
  R.version$major,
  strsplit(R.version$minor, "\\.")[[1]][1],
  sep = "."
)
project_library <- file.path(script_dir, "R_library", r_major_minor)
legacy_library <- file.path(script_dir, "R-library")
shared_project_library <- file.path(dirname(script_dir), "R-library")
dir.create(project_library, recursive = TRUE, showWarnings = FALSE)
.libPaths(unique(c(
  project_library,
  shared_project_library[file.exists(shared_project_library)],
  legacy_library[file.exists(legacy_library)],
  .libPaths()
)))
if (identical(getOption("repos")["CRAN"], "@CRAN@") ||
    is.na(getOption("repos")["CRAN"])) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}
required_packages <- c(
  "dplyr", "tidyr", "purrr", "tibble", "ggplot2", "ggrepel",
  "metafor", "clubSandwich", "svglite", "ragg", "scales", "patchwork"
)
optional_packages <- c(
  "brms", "posterior"
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
install_or_stop(required_packages, required = TRUE, auto_install = FALSE)
missing_optional <- install_or_stop(optional_packages, required = FALSE, auto_install = FALSE)
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(metafor)
  library(clubSandwich)
  library(scales)
  library(patchwork)
})
# ------------------------------ 1. Figure contract ----------------------------
# Core conclusion:
# Unfavourable tumor-informed ctDNA status is associated with worse time-to-event
# outcomes, with prognostic strength depending on sampling timepoint.
#
# Evidence chain:
# 1) timepoint-specific recurrence-related outcomes;
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
  if (!isTRUE(MAKE_FIGURES)) return(invisible(NULL))
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
  if (!isTRUE(MAKE_FIGURES)) return(invisible(NULL))
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
  if (!isTRUE(MAKE_FIGURES)) return(invisible(NULL))
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
# Embedded source data copied verbatim from the final extraction workbook.
# Effect_Estimates: 92 candidate effect estimates used to construct all
# prespecified primary, adjusted, OS, subtype, and sensitivity datasets.
raw_effects <- tibble::as_tibble(
structure(list(Effect_ID = c("E001", "E002", "E003", "E004", 
"E005", "E006", "E007", "E008", "E009", "E010", "E011", "E083", 
"E084", "E012", "E013", "E014", "E015", "E016", "E017", "E018", 
"E019", "E020", "E021", "E022", "E085", "E086", "E023", "E024", 
"E025", "E026", "E079", "E080", "E081", "E082", "E027", "E028", 
"E029", "E030", "E031", "E032", "E033", "E034", "E035", "E036", 
"E037", "E038", "E039", "E040", "E041", "E042", "E043", "E044", 
"E045", "E046", "E047", "E048", "E049", "E050", "E051", "E052", 
"E053", "E054", "E055", "E056", "E057", "E058", "E059", "E060", 
"E061", "E062", "E063", "E064", "E065", "E066", "E067", "E068", 
"E069", "E070", "E071", "E072", "E073", "E087", "E088", "E089", 
"E090", "E091", "E092", "E074", "E075", "E076", "E077", "E078"
), Study_ID = c("S01", "S01", "S01", "S01", "S02", "S03", "S03", 
"S03", "S03", "S03", "S03", "S03", "S03", "S04", "S04", "S04", 
"S04", "S04", "S04", "S04", "S04", "S04", "S04", "S04", "S04", 
"S04", "S05", "S05", "S05", "S05", "S05", "S05", "S05", "S05", 
"S06", "S07", "S07", "S07", "S07", "S07", "S07", "S07", "S07", 
"S08", "S08", "S08", "S08", "S08", "S08", "S08", "S08", "S09", 
"S09", "S09", "S09", "S09", "S09", "S10", "S10", "S10", "S10", 
"S10", "S11", "S11", "S11", "S11", "S11", "S11", "S12", "S12", 
"S12", "S12", "S12", "S12", "S12", "S12", "S12", "S13", "S13", 
"S13", "S13", "S13", "S13", "S13", "S13", "S13", "S13", "S14", 
"S14", "S14", "S14", "S14"), Citation = c("Cavallone 2020", "Cavallone 2020", 
"Cavallone 2020", "Cavallone 2020", "Garcia-Murillas 2025", "Elliott 2025", 
"Elliott 2025", "Elliott 2025", "Elliott 2025", "Elliott 2025", 
"Elliott 2025", "Elliott 2025", "Elliott 2025", "Li 2025", "Li 2025", 
"Li 2025", "Li 2025", "Li 2025", "Li 2025", "Li 2025", "Li 2025", 
"Li 2025", "Li 2025", "Li 2025", "Li 2025", "Li 2025", "IMpassion031 2025", 
"IMpassion031 2025", "IMpassion031 2025", "IMpassion031 2025", 
"IMpassion031 2025", "IMpassion031 2025", "IMpassion031 2025", 
"IMpassion031 2025", "Cabel 2026 (RaDaR)", "Cabel 2026 (SCANDARE)", 
"Cabel 2026 (SCANDARE)", "Cabel 2026 (SCANDARE)", "Cabel 2026 (SCANDARE)", 
"Cabel 2026 (SCANDARE)", "Cabel 2026 (SCANDARE)", "Cabel 2026 (SCANDARE)", 
"Cabel 2026 (SCANDARE)", "George 2026", "George 2026", "George 2026", 
"George 2026", "George 2026", "George 2026", "George 2026", "George 2026", 
"Blaye 2026", "Blaye 2026", "Blaye 2026", "Blaye 2026", "Blaye 2026", 
"Blaye 2026", "Roseshter 2026", "Roseshter 2026", "Roseshter 2026", 
"Roseshter 2026", "Roseshter 2026", "Hunter 2026", "Hunter 2026", 
"Hunter 2026", "Hunter 2026", "Hunter 2026", "Hunter 2026", "Egle 2026", 
"Egle 2026", "Egle 2026", "Egle 2026", "Egle 2026", "Egle 2026", 
"Egle 2026", "Egle 2026", "Egle 2026", "Magbanua 2025", "Magbanua 2025", 
"Magbanua 2025", "Magbanua 2025", "Magbanua 2025", "Magbanua 2025", 
"Magbanua 2025", "Magbanua 2025", "Magbanua 2025", "Magbanua 2025", 
"Cailleux 2022", "Cailleux 2022", "Cailleux 2022", "Cailleux 2022", 
"Cailleux 2022"), Cohort = c("Q-CROC-03", "Q-CROC-03", "Q-CROC-03", 
"Q-CROC-03", "ChemoNEAR", "LIBERATE", "LIBERATE", "LIBERATE", 
"LIBERATE", "LIBERATE", "LIBERATE", "LIBERATE", "LIBERATE", "Primary TNBC cohort", 
"Primary TNBC cohort", "Primary TNBC cohort", "Primary TNBC cohort", 
"Primary TNBC cohort", "Primary TNBC cohort", "Primary TNBC cohort", 
"I-SPY2 external validation", "Primary TNBC cohort", "Primary TNBC cohort", 
"Primary TNBC cohort", "Primary TNBC cohort", "Primary TNBC cohort", 
"IMpassion031", "IMpassion031", "IMpassion031", "IMpassion031", 
"IMpassion031", "IMpassion031", "IMpassion031", "IMpassion031", 
"MSK cohort", "SCANDARE", "SCANDARE", "SCANDARE", "SCANDARE", 
"SCANDARE", "SCANDARE", "SCANDARE", "SCANDARE", "NeoCircle", 
"NeoCircle", "NeoCircle", "NeoCircle", "NeoCircle", "NeoCircle", 
"NeoCircle", "NeoCircle", "ALIENOR", "ALIENOR", "ALIENOR", "ALIENOR", 
"ALIENOR", "ALIENOR", "TRICIA/JGH", "TRICIA/JGH", "TRICIA/JGH", 
"TRICIA/JGH", "TRICIA/JGH", "PREDICT-DNA", "PREDICT-DNA", "PREDICT-DNA", 
"PREDICT-DNA", "PREDICT-DNA", "PREDICT-DNA", "ABCSG-34", "ABCSG-34", 
"ABCSG-34", "ABCSG-34", "ABCSG-34", "ABCSG-34", "ABCSG-34", "ABCSG-34", 
"ABCSG-34", "I-SPY2", "I-SPY2", "I-SPY2", "I-SPY2", "I-SPY2", 
"I-SPY2", "I-SPY2", "I-SPY2", "I-SPY2", "I-SPY2", "Single-center cohort", 
"Single-center cohort", "Single-center cohort", "Single-center cohort", 
"Single-center cohort"), Timepoint = c("After cycle 1", "After cycle 1", 
"End-NAC/preop", "End-NAC/preop", "Any postsurgery surveillance", 
"Baseline", "Baseline", "Pre-cycle 2", "Mid-NAT", "Mid-NAT", 
"Any postoperative/follow-up", "Mid-NAT", "Mid-NAT", "Baseline", 
"Post-NAC/preop", "Post-NAC/preop", "Postsurgery", "Postsurgery", 
"Baseline MVAF >1.1%", "Baseline MVAF >1.1%", "Baseline high MTM/mL", 
"Composite systemic burden", "Composite systemic burden", "Longitudinal surveillance", 
"Composite systemic burden", "Longitudinal surveillance", "Week 7 during NAT", 
"Week 7 during NAT", "Week 7 during NAT", "Week 7 during NAT", 
"Longitudinal through week 7", "Longitudinal through week 7", 
"Post-surgery", "Post-surgery", "Any postsurgery/follow-up", 
"Post-NAT/preop", "Post-NAT/preop", "Post-NAT/preop", "Post-NAT/preop continuous", 
"Baseline continuous", "During-NAT dynamics", "During-NAT dynamics", 
"During-NAT dynamics", "End-NAT/preop", "End-NAT/preop", "NAT dynamics", 
"NAT dynamics", "Postoperative landmark", "Postoperative landmark", 
"Any follow-up MRD", "Any follow-up MRD", "Any postsurgery follow-up", 
"Any postsurgery follow-up", "First postsurgery sample", "First postsurgery sample", 
"First postsurgery sample", "First postsurgery sample", "Post-NAC/preop T1", 
"Post-NAC/preop T1", "Post-NAC/preop T1", "Late postsurgery T4", 
"Late postsurgery T4", "Post-NAT/preop T1", "Post-NAT/preop T1", 
"Post-NAT/preop T1", "Post-NAT/preop T1", "Postsurgery T2*", 
"Postsurgery T2*", "Baseline", "Baseline", "Baseline", "Mid-therapy", 
"Mid-therapy", "Mid-therapy", "End-of-treatment", "End-of-treatment", 
"End-of-treatment", "Pretreatment T0", "Post-NAT/preop T3", "Clearance trajectory", 
"Clearance trajectory", "Clearance trajectory", "Clearance trajectory", 
"Pretreatment T0", "Pretreatment T0", "Post-NAT/preop T3", "Post-NAT/preop T3", 
"Baseline", "Presurgery", "Presurgery", "Last follow-up", "Last follow-up"
), `Analysis N` = c(26, 26, 26, 26, 61, 114, 41, 82, 74, 45, 
86, 26, 19, 122, 113, 113, 118, 118, 122, 122, NA, 119, 119, 
82, 119, 82, 121, 123, 68, 68, 130, 130, 120, 120, 30, 68, 68, 
NA, 68, 86, NA, NA, NA, 131, 131, 122, 122, 134, 134, 136, 136, 
83, 83, 83, 83, 83, 83, 64, 64, 64, 66, 66, 129, 68, 68, 61, 
75, 76, 109, 109, 109, 95, 95, 95, 89, 89, 89, 712, 712, NA, 
NA, NA, NA, NA, NA, NA, NA, 38, 41, 41, 38, 38), Subgroup = c("All", 
"All", "All", "All", "All", "All", "ER+", "All", "All", "HER2-negative", 
"All", "ER+", "TNBC", "All", "All", "All", "All", "All", "All", 
"All", "All", "All", "All", "All", "All", "All", "All", "All", 
"Non-pCR", "Non-pCR", "All", "All", "All", "All", "All", "All", 
"Ultrasensitive <100 ppm", "Non-pCR", "All", "All", "All", "All", 
"All", "All", "All", "All", "All", "All", "All", "All", "All", 
"All", "All", "All", "All", "All", "All", "TNBC non-pCR", "TNBC non-pCR", 
"TNBC non-pCR", "TNBC non-pCR", "TNBC non-pCR", "TNBC + HER2+", 
"TNBC", "TNBC", "HER2+", "TNBC", "HER2+", "All", "All", "All", 
"All", "All", "All", "All", "All", "All", "All", "All", "All", 
"All", "All", "All", "RCB-II", "RCB-III", "RCB-II", "RCB-III", 
"All", "All", "All", "All", "All"), Endpoint = c("RFS", "OS", 
"RFS", "OS", "RFS", "RFI", "RFI", "RFI", "RFI", "RFI", "RFI", 
"RFI", "RFI", "EFS", "EFS", "DRFS", "EFS", "DRFS", "EFS", "DRFS", 
"DRFS", "EFS", "DRFS", "EFS", "DRFS", "EFS", "DFS", "OS", "DFS", 
"OS", "DFS", "OS", "DFS", "OS", "DFS", "DRFI", "DRFI", "DRFI", 
"DRFI", "DRFI", "DRFI", "DFS", "OS", "BCFi", "OS", "BCFi", "OS", 
"BCFi", "OS", "BCFi", "OS", "RFI", "RFI", "RFI", "RFI", "OS", 
"OS", "RFS", "OS", "RFS", "RFS", "OS", "IDFS", "IDFS", "IDFS", 
"IDFS", "IDFS", "IDFS", "iDFS", "DRFS", "OS", "iDFS", "DRFS", 
"OS", "iDFS", "DRFS", "OS", "DRFS", "DRFS", "DRFS", "DRFS", "DRFS", 
"DRFS", "DRFS", "DRFS", "DRFS", "DRFS", "EFS", "EFS", "EFS", 
"EFS", "EFS"), Analysis = c("Univariable", "Univariable", "Univariable", 
"Univariable", "Time-dependent univariable Cox", "Univariable", 
"Univariable", "Univariable", "Univariable", "Univariable", "Univariable", 
"Univariable", "Univariable", "Univariable", "Univariable", "Univariable", 
"Univariable", "Univariable", "Univariable", "Univariable", "Univariable", 
"Univariable", "Univariable", "Univariable", "Multivariable", 
"Multivariable", "Univariable", "Univariable", "Univariable", 
"Univariable", "Univariable", "Univariable", "Univariable", "Univariable", 
"Univariable", "Univariable", "Univariable", "Univariable", "Multivariable", 
"Multivariable", "Univariable", "Univariable", "Univariable", 
"Univariable", "Univariable", "Univariable", "Univariable", "Univariable", 
"Univariable", "Univariable", "Univariable", "Univariable", "Multivariable", 
"Univariable", "Multivariable", "Univariable", "Multivariable", 
"Univariable", "Univariable", "Multivariable", "Univariable", 
"Univariable", "Multivariable", "Univariable", "Multivariable", 
"Univariable", "Univariable", "Univariable", "Univariable", "Univariable", 
"Univariable", "Univariable", "Univariable", "Univariable", "Univariable", 
"Univariable", "Univariable", "Multivariable", "Multivariable", 
"Multivariable", "Multivariable", "Multivariable", "Multivariable", 
"Adjusted/stratified Cox", "Adjusted/stratified Cox", "Adjusted/stratified Cox", 
"Adjusted/stratified Cox", "Univariable", "Univariable", "Multivariable", 
"Univariable", "Multivariable"), `Reported contrast` = c("ctDNA-negative vs ctDNA-positive", 
"ctDNA-negative vs ctDNA-positive", "ctDNA-negative vs ctDNA-positive", 
"ctDNA-negative vs ctDNA-positive", "Any ctDNA-positive vs always negative", 
"ctDNA-positive vs negative", "ctDNA-positive vs negative", "Persistent ctDNA-positive vs cleared/negative", 
"ctDNA-positive vs negative", "ctDNA-positive vs negative", "Any ctDNA-positive vs always negative", 
"ctDNA-positive vs negative", "ctDNA-positive vs negative", "ctDNA-positive vs negative", 
"ctDNA-positive vs negative", "ctDNA-positive vs negative", "ctDNA-positive vs negative", 
"ctDNA-positive vs negative", "High vs low MVAF", "High vs low MVAF", 
"High vs low ctDNA burden", "High vs low composite burden", "High vs low composite burden", 
"MRD-positive vs MRD-negative", "High vs low composite burden", 
"MRD-positive vs MRD-negative", "Clearance vs no clearance", 
"Clearance vs no clearance", "ctDNA-positive vs negative", "ctDNA-positive vs negative", 
"Always ctDNA-positive vs negative at <U+2265>1 timepoint", "Always ctDNA-positive vs negative at <U+2265>1 timepoint", 
"ctDNA-positive vs negative", "ctDNA-positive vs negative", "Any ctDNA-positive vs always negative", 
"ctDNA-positive vs negative", "Ultrasensitive positive vs negative", 
"ctDNA-positive vs negative", "Per 10-fold ctDNA increase", "Per ctDNA level unit reported", 
"Increase vs decrease", "Increase vs decrease", "Increase vs decrease", 
"ctDNA-positive vs negative", "ctDNA-positive vs negative", "Nonresponse vs response", 
"Nonresponse vs response", "ctDNA-positive vs negative", "ctDNA-positive vs negative", 
"Any ctDNA-positive vs always negative", "Any ctDNA-positive vs always negative", 
"Any ctDNA-positive vs always negative", "Any ctDNA-positive vs always negative", 
"ctDNA-positive vs negative", "ctDNA-positive vs negative", "ctDNA-positive vs negative", 
"ctDNA-positive vs negative", "ctDNA-negative vs positive", "ctDNA-negative vs positive", 
"ctDNA-positive vs negative", "ctDNA-negative vs positive", "ctDNA-negative vs positive", 
"ctDNA-positive vs negative", "ctDNA-positive vs negative", "ctDNA-positive vs negative", 
"ctDNA-positive vs negative", "ctDNA-positive vs negative", "ctDNA-positive vs negative", 
"ctDNA-positive vs negative", "ctDNA-positive vs negative", "ctDNA-positive vs negative", 
"ctDNA-positive vs negative", "ctDNA-positive vs negative", "ctDNA-positive vs negative", 
"ctDNA-positive vs negative", "ctDNA-positive vs negative", "ctDNA-positive vs negative", 
"ctDNA-positive vs negative", "ctDNA-positive vs negative", "Cleared at T3 vs persistently negative", 
"No clearance vs persistently negative", "Cleared at T1 (week 3) vs persistently negative", 
"Cleared at T2 (week 12) vs persistently negative", "ctDNA-positive vs negative", 
"ctDNA-positive vs negative", "ctDNA-positive vs negative", "ctDNA-positive vs negative", 
"ctDNA-positive vs negative", "ctDNA-positive vs negative", "ctDNA-positive vs negative", 
"ctDNA-positive vs negative", "ctDNA-positive vs negative"), 
    HR_reported = c(0.32, 0.35, 0.29, 0.27, 37.2, 2.89, 3.03, 
    1.13, 3.59, 14.43, 35.7, 7.59, 30.1, 3.39, 4.34, 6.18, 4.63, 
    6.25, 4.49, 13.02, 2.55, 3.76, 11.04, 7.69, 11.51, 7.99, 
    0.26, 0.2, 2.15, 3.5, 3.6, 4.57, 40, 107.74, 52, 35.6, 32.1, 
    14.8, 2, 1.7, 7.3, 8.2, 5.5, 3.3, 5.6, 9.4, 11.9, 3.5, 5.6, 
    37.9, 12.7, 29.2, 69.7, 5.3, 4.6, 3.9, 5.9, 0.052, 0.06, 
    24.7, 0.24, 0.147, 8.9, 8.9, 12, 4.4, 128, 263, 1.53, 1.98, 
    2.12, 1.35, 1.71, 1.9, 1.09, 1.09, 1.27, 4.4, 5.2, 6.88, 
    16.5, 2.03, 1.94, 10.1, 29.97, 3.47, 13.22, 1.4, 61, 53, 
    23, 31), CI_low_reported = c(0.09, 0.09, 0.08, 0.075, 10.5, 
    1.003, 0.87, 0.44, 0.99, 3, 4.3, 0.99, 2.65, 0.44, 1.51, 
    1.95, 1.6, 2.01, 1.21, 1.65, 1.37, 1.05, 1.42, 1.7, 1.48, 
    1.73, 0.12, 0.07, 0.88, 1.21, 1.36, 1.48, 10.41, 19.1, 5, 
    4.6, 3.8, 1.9, 1.4, 0.9, 0.9, 1, 1.2, 1.4, 1.9, 3.8, 3.9, 
    1.2, 1.7, 12.7, 4.1, 6.9, 13.3, 1.1, 1.8, 1.2, 1.6, 0.024, 
    0.027, 3.2, 0.093, 0.037, 1.92, 2.4, 1.52, 1, 15, 27.2, 0.81, 
    0.99, 1.02, 0.62, 0.76, 0.84, 0.44, 0.41, 0.47, 1.91, 3.24, 
    2.37, 5.67, 0.67, 0.6, 1.38, 4.09, 1.63, 7.47, 0.32, 5.4, 
    4.5, 3.1, 2.7), CI_high_reported = c(1.1, 1.36, 0.98, 0.96, 
    131.9, 8.31, 10.62, 2.96, 12.87, 69.53, 297, 57.86, 342.4, 
    26.3, 12.44, 19.58, 13.38, 19.42, 16.6, 102.8, 4.73, 13.49, 
    85.5, 34.74, 89.22, 36.85, 0.57, 0.54, 5.29, 10.16, 9.57, 
    14.05, 153.69, 607.66, 522, 275, 268.2, 114.3, 2.8, 3.1, 
    58.3, 64.5, 25.4, 7.9, 16.7, 23.7, 35.9, 10.6, 18.5, 113.6, 
    39.2, 124, 365.1, 5.7, 12.1, 12.8, 22.4, 0.109, 0.131, 191.34, 
    0.62, 0.59, 41.1, 33, 95, 18, 1083, 3538, 2.86, 3.95, 4.41, 
    2.93, 3.83, 4.29, 2.69, 2.91, 3.44, 10.16, 8.35, 19.98, 47.98, 
    6.19, 6.25, 74.09, 219.69, 7.41, 23.37, 5.9, 689, 624, 170, 
    352), `Invert?` = c("Yes", "Yes", "Yes", "Yes", "No", "No", 
    "No", "No", "No", "No", "No", "No", "No", "No", "No", "No", 
    "No", "No", "No", "No", "No", "No", "No", "No", "No", "No", 
    "Yes", "Yes", "No", "No", "No", "No", "No", "No", "No", "No", 
    "No", "No", "No", "No", "No", "No", "No", "No", "No", "No", 
    "No", "No", "No", "No", "No", "No", "No", "No", "No", "No", 
    "No", "Yes", "Yes", "No", "Yes", "Yes", "No", "No", "No", 
    "No", "No", "No", "No", "No", "No", "No", "No", "No", "No", 
    "No", "No", "No", "No", "No", "No", "No", "No", "No", "No", 
    "No", "No", "No", "No", "No", "No", "No"), HR_standardized = c(3.125, 
    2.85714285714286, 3.44827586206897, 3.7037037037037, 37.2, 
    2.89, 3.03, 1.13, 3.59, 14.43, 35.7, 7.59, 30.1, 3.39, 4.34, 
    6.18, 4.63, 6.25, 4.49, 13.02, 2.55, 3.76, 11.04, 7.69, 11.51, 
    7.99, 3.84615384615385, 5, 2.15, 3.5, 3.6, 4.57, 40, 107.74, 
    52, 35.6, 32.1, 14.8, 2, 1.7, 7.3, 8.2, 5.5, 3.3, 5.6, 9.4, 
    11.9, 3.5, 5.6, 37.9, 12.7, 29.2, 69.7, 5.3, 4.6, 3.9, 5.9, 
    19.2307692307692, 16.6666666666667, 24.7, 4.16666666666667, 
    6.80272108843537, 8.9, 8.9, 12, 4.4, 128, 263, 1.53, 1.98, 
    2.12, 1.35, 1.71, 1.9, 1.09, 1.09, 1.27, 4.4, 5.2, 6.88, 
    16.5, 2.03, 1.94, 10.1, 29.97, 3.47, 13.22, 1.4, 61, 53, 
    23, 31), CI_low_standardized = c(0.909090909090909, 0.735294117647059, 
    1.02040816326531, 1.04166666666667, 10.5, 1.003, 0.87, 0.44, 
    0.99, 3, 4.3, 0.99, 2.65, 0.44, 1.51, 1.95, 1.6, 2.01, 1.21, 
    1.65, 1.37, 1.05, 1.42, 1.7, 1.48, 1.73, 1.75438596491228, 
    1.85185185185185, 0.88, 1.21, 1.36, 1.48, 10.41, 19.1, 5, 
    4.6, 3.8, 1.9, 1.4, 0.9, 0.9, 1, 1.2, 1.4, 1.9, 3.8, 3.9, 
    1.2, 1.7, 12.7, 4.1, 6.9, 13.3, 1.1, 1.8, 1.2, 1.6, 9.1743119266055, 
    7.63358778625954, 3.2, 1.61290322580645, 1.69491525423729, 
    1.92, 2.4, 1.52, 1, 15, 27.2, 0.81, 0.99, 1.02, 0.62, 0.76, 
    0.84, 0.44, 0.41, 0.47, 1.91, 3.24, 2.37, 5.67, 0.67, 0.6, 
    1.38, 4.09, 1.63, 7.47, 0.32, 5.4, 4.5, 3.1, 2.7), CI_high_standardized = c(11.1111111111111, 
    11.1111111111111, 12.5, 13.3333333333333, 131.9, 8.31, 10.62, 
    2.96, 12.87, 69.53, 297, 57.86, 342.4, 26.3, 12.44, 19.58, 
    13.38, 19.42, 16.6, 102.8, 4.73, 13.49, 85.5, 34.74, 89.22, 
    36.85, 8.33333333333333, 14.2857142857143, 5.29, 10.16, 9.57, 
    14.05, 153.69, 607.66, 522, 275, 268.2, 114.3, 2.8, 3.1, 
    58.3, 64.5, 25.4, 7.9, 16.7, 23.7, 35.9, 10.6, 18.5, 113.6, 
    39.2, 124, 365.1, 5.7, 12.1, 12.8, 22.4, 41.6666666666667, 
    37.037037037037, 191.34, 10.752688172043, 27.027027027027, 
    41.1, 33, 95, 18, 1083, 3538, 2.86, 3.95, 4.41, 2.93, 3.83, 
    4.29, 2.69, 2.91, 3.44, 10.16, 8.35, 19.98, 47.98, 6.19, 
    6.25, 74.09, 219.69, 7.41, 23.37, 5.9, 689, 624, 170, 352
    ), logHR = c(1.13943428318836, 1.04982212449868, 1.23787435600162, 
    1.30933331998376, 3.6163087612791, 1.06125650212434, 1.10856261952128, 
    0.122217632724249, 1.27815220250019, 2.66930937278578, 3.57515068878559, 
    2.02683159140754, 3.40452517175483, 1.22082992139236, 1.46787434811231, 
    1.8213182714696, 1.53255686809814, 1.83258146374831, 1.50185270175416, 
    2.56648663678042, 0.936093359170335, 1.3244189574018, 2.40152504084895, 
    2.03992078351755, 2.44321622273379, 2.07819075977818, 1.34707364796661, 
    1.6094379124341, 0.765467842139571, 1.25276296849537, 1.28093384546206, 
    1.51951320490611, 3.68887945411394, 4.67972091725239, 3.95124371858143, 
    3.57234563785798, 3.46885603013597, 2.69462718077007, 0.693147180559945, 
    0.53062825106217, 1.98787434815435, 2.10413415427021, 1.70474809223843, 
    1.19392246847243, 1.7227665977411, 2.24070968927596, 2.47653840011748, 
    1.25276296849537, 1.7227665977411, 3.63495111208838, 2.54160199346455, 
    3.37416870927424, 4.24420031776648, 1.66770682055808, 1.52605630349505, 
    1.3609765531356, 1.77495235091167, 2.95651156040071, 2.81341071676004, 
    3.20680324363393, 1.42711635564015, 1.9173226922034, 2.18605127673809, 
    2.18605127673809, 2.484906649788, 1.48160454092422, 4.85203026391962, 
    5.57215403217776, 0.425267735404344, 0.683096844706444, 0.751416088683921, 
    0.300104592450338, 0.536493370514568, 0.641853886172395, 
    0.0861776962410524, 0.0861776962410524, 0.2390169004705, 
    1.48160454092422, 1.64865862558738, 1.92861865194525, 2.80336038090653, 
    0.708035793053696, 0.662687973075237, 2.31253542384721, 3.40019688132857, 
    1.24415459395877, 2.58173083442354, 0.336472236621213, 4.11087386417331, 
    3.97029191355212, 3.13549421592915, 3.43398720448515), SE_logHR = c(0.638585660320458, 
    0.692711813367304, 0.639164779844576, 0.650368666052442, 
    0.645578776163212, 0.539404107114173, 0.638265582435536, 
    0.48626781132801, 0.654323815678963, 0.801822965674188, 1.08038701941408, 
    1.03777465188814, 1.24015816435445, 1.04350752327897, 0.537961080735765, 
    0.588438587889467, 0.541774853432231, 0.578614474691886, 
    0.668056718304553, 1.0540842002828, 0.316100628230563, 0.651315945420978, 
    1.04537232253305, 0.769710079018058, 1.04567937293982, 0.78028933931777, 
    0.39748587195065, 0.521192320793102, 0.457564188107167, 0.542815837383083, 
    0.497741965744108, 0.574127629592789, 0.686778240814373, 
    0.882634484686414, 1.18577287638993, 1.04354969239069, 1.08590099594776, 
    1.04514609301744, 0.176823260346925, 0.315500670191053, 1.06402107372813, 
    1.06292480199024, 0.778680514599627, 0.441426153788715, 0.554478273788346, 
    0.466957648318558, 0.566265495499651, 0.55574807253165, 0.608964918628089, 
    0.558949365515945, 0.575941268644427, 0.736928610714751, 
    0.845001850365786, 0.419682651794944, 0.486076221352188, 
    0.603858064829494, 0.673228910616137, 0.386044400990102, 
    0.402898075815516, 1.04359722970519, 0.483959179817827, 0.706429751085597, 
    0.781551769250638, 0.668632353089944, 1.0548894277404, 0.737339734157185, 
    1.09169388925989, 1.24186221974943, 0.321822106160055, 0.353001508869013, 
    0.373487770959841, 0.396183220401014, 0.41257695124855, 0.415979622470588, 
    0.461867792266219, 0.499936530731425, 0.507779095832485, 
    0.426366122472397, 0.241502094147583, 0.543837198488104, 
    0.544794678487726, 0.567197105432048, 0.597807930488342, 
    1.0161217018519, 1.01624298257306, 0.386288373594215, 0.290958711583174, 
    0.743465978086745, 1.23694957077639, 1.25818188050846, 1.02152967488754, 
    1.24244372514995), Variance = c(0.407791645566916, 0.479849656378619, 
    0.408531615793765, 0.422979401782832, 0.416771956232391, 
    0.290956790771638, 0.407382953721774, 0.236456384333733, 
    0.428139655764678, 0.64292006828255, 1.16723611171843, 1.07697622810154, 
    1.53799227261501, 1.0889079511398, 0.289402124386393, 0.346259971717349, 
    0.293519991811516, 0.334794710322968, 0.446299778871849, 
    1.11109350128583, 0.0999196071677568, 0.424212460759622, 
    1.09280329271814, 0.592453605741986, 1.09344535099182, 0.608851453052962, 
    0.157995018400369, 0.2716414352537, 0.209364986238171, 0.294649033313898, 
    0.247747064462809, 0.329622535061834, 0.471664352056085, 
    0.779043633557652, 1.40605731438206, 1.08899596048871, 1.17918097300035, 
    1.09233035574962, 0.0312664653997164, 0.0995406728910034, 
    1.13214084533756, 1.12980913468598, 0.606343343817139, 0.194857049248698, 
    0.307446156103304, 0.218049445323198, 0.320656611393465, 
    0.308855920122644, 0.370838272119715, 0.312424393210678, 
    0.331708344927752, 0.543063777289973, 0.714028127121603, 
    0.176133528217636, 0.236270092964021, 0.364644562459622, 
    0.453237166089391, 0.149030279535806, 0.162326859495845, 
    1.08909517784834, 0.234216487729944, 0.499042993218859, 0.610823168018802, 
    0.447069223598595, 1.11279170475846, 0.543669883566988, 1.19179554784739, 
    1.54222177284099, 0.103569468013293, 0.1246100652638, 0.13949311505655, 
    0.156961144127319, 0.170219740701548, 0.173039046310773, 
    0.213321857532871, 0.249936534759773, 0.257839610164456, 
    0.181788070392147, 0.0583232614776681, 0.29575889845939, 
    0.296801241708545, 0.321712556410494, 0.357374321754754, 
    1.0325033129744, 1.03274979962898, 0.149218707574064, 0.0846569718461404, 
    0.55274166057248, 1.53004424064388, 1.5830216444398, 1.04352287667585, 
    1.54366641016449), `Recommended use` = c("No<U+2014>exploratory early timepoint", 
    "No<U+2014>exploratory early timepoint", "Yes<U+2014>preoperative RFS model", 
    "Secondary<U+2014>OS", "Yes<U+2014>postoperative surveillance", 
    "Yes<U+2014>baseline RFI", "Subtype only", "Early-treatment subgroup", 
    "Yes<U+2014>mid-NAT, sensitivity", "Subtype only", "Yes<U+2014>postoperative surveillance", 
    "Subtype only", "Subtype only", "Baseline sensitivity only", 
    "Yes<U+2014>preoperative EFS", "Yes<U+2014>preoperative DRFS", 
    "Yes<U+2014>postoperative EFS", "Yes<U+2014>postoperative DRFS", 
    "Continuous/burden subgroup only", "Continuous/burden subgroup only", 
    "Do not pool<U+2014>overlaps S13", "Composite predictor only", 
    "Composite predictor only", "Surveillance sensitivity", "Adjusted sensitivity", 
    "Adjusted sensitivity", "Yes<U+2014>week-7 clearance DFS", 
    "Secondary<U+2014>week-7 OS", "Non-pCR subgroup only", "Non-pCR subgroup only", 
    "Yes<U+2014>longitudinal persistence DFS", "Secondary<U+2014>longitudinal OS", 
    "Yes<U+2014>post-surgery DFS", "Secondary<U+2014>post-surgery OS", 
    "Yes<U+2014>postoperative surveillance", "Yes<U+2014>preoperative DRFI", 
    "Assay-sensitivity subgroup", "Non-pCR subgroup", "Continuous effect only", 
    "Continuous effect only", "Dynamics subgroup", "Dynamics subgroup", 
    "Dynamics subgroup", "Yes<U+2014>preoperative BCFi", "Secondary<U+2014>OS", 
    "Yes<U+2014>dynamics BCFi", "Secondary<U+2014>dynamics OS", 
    "Yes<U+2014>postoperative landmark", "Secondary<U+2014>OS", 
    "Yes<U+2014>surveillance BCFi", "Secondary<U+2014>OS", "Yes<U+2014>surveillance RFI", 
    "Adjusted sensitivity", "Hold<U+2014>verify CI", "Yes<U+2014>postoperative adjusted", 
    "Secondary<U+2014>OS", "Secondary<U+2014>adjusted OS", "Yes<U+2014>preoperative RFS", 
    "Secondary<U+2014>OS", "Adjusted sensitivity", "Yes<U+2014>late postoperative RFS", 
    "Secondary<U+2014>OS", "Yes<U+2014>adjusted preoperative IDFS", 
    "Subtype primary", "Subtype adjusted", "Subtype exploratory", 
    "Postoperative subgroup; sparse", "Postoperative subgroup; sparse", 
    "Yes<U+2014>baseline iDFS", "Yes<U+2014>baseline DRFS", "Secondary<U+2014>OS", 
    "Yes<U+2014>mid-NAT iDFS", "Yes<U+2014>mid-NAT DRFS", "Secondary<U+2014>OS", 
    "Yes<U+2014>preoperative iDFS", "Yes<U+2014>preoperative DRFS", 
    "Secondary<U+2014>OS", "Yes<U+2014>adjusted baseline DRFS", 
    "Yes<U+2014>adjusted preoperative DRFS", "Dynamics subgroup", 
    "Dynamics subgroup", "Dynamics subgroup", "Dynamics subgroup", 
    "RCB subgroup only", "RCB subgroup only", "RCB subgroup only", 
    "RCB subgroup only", "Yes<U+2014>baseline EFS", "Yes<U+2014>preoperative EFS", 
    "Adjusted sensitivity", "Yes<U+2014>surveillance EFS", "Adjusted sensitivity"
    ), `Dependency cluster` = c("S01_QCROC", "S01_QCROC", "S01_QCROC", 
    "S01_QCROC", "S02_ChemoNEAR", "S03_LIBERATE", "S03_LIBERATE", 
    "S03_LIBERATE", "S03_LIBERATE", "S03_LIBERATE", "S03_LIBERATE", 
    "S03_LIBERATE", "S03_LIBERATE", "S04_LiPrimary", "S04_LiPrimary", 
    "S04_LiPrimary", "S04_LiPrimary", "S04_LiPrimary", "S04_LiPrimary", 
    "S04_LiPrimary", "S13_ISPY2", "S04_LiPrimary", "S04_LiPrimary", 
    "S04_LiPrimary", "S04_LiPrimary", "S04_LiPrimary", "S05_IMpassion031", 
    "S05_IMpassion031", "S05_IMpassion031", "S05_IMpassion031", 
    "S05_IMpassion031", "S05_IMpassion031", "S05_IMpassion031", 
    "S05_IMpassion031", "S06_MSK_RaDaR", "S07_SCANDARE", "S07_SCANDARE", 
    "S07_SCANDARE", "S07_SCANDARE", "S07_SCANDARE", "S07_SCANDARE", 
    "S07_SCANDARE", "S07_SCANDARE", "S08_NeoCircle", "S08_NeoCircle", 
    "S08_NeoCircle", "S08_NeoCircle", "S08_NeoCircle", "S08_NeoCircle", 
    "S08_NeoCircle", "S08_NeoCircle", "S09_ALIENOR", "S09_ALIENOR", 
    "S09_ALIENOR", "S09_ALIENOR", "S09_ALIENOR", "S09_ALIENOR", 
    "S10_TRICIA", "S10_TRICIA", "S10_TRICIA", "S10_TRICIA", "S10_TRICIA", 
    "S11_PREDICT", "S11_PREDICT", "S11_PREDICT", "S11_PREDICT", 
    "S11_PREDICT", "S11_PREDICT", "S12_ABCSG34", "S12_ABCSG34", 
    "S12_ABCSG34", "S12_ABCSG34", "S12_ABCSG34", "S12_ABCSG34", 
    "S12_ABCSG34", "S12_ABCSG34", "S12_ABCSG34", "S13_ISPY2", 
    "S13_ISPY2", "S13_ISPY2", "S13_ISPY2", "S13_ISPY2", "S13_ISPY2", 
    "S13_ISPY2", "S13_ISPY2", "S13_ISPY2", "S13_ISPY2", "S14_Cailleux", 
    "S14_Cailleux", "S14_Cailleux", "S14_Cailleux", "S14_Cailleux"
    ), `Source locator` = c("Results/Fig. 5", "Results/Fig. 5", 
    "Abstract/Results", "Abstract/Results", "Results/Fig. 3", 
    "Results/Fig. 4", "Supplementary Fig. 3A", "Supplementary Fig. 3C", 
    "Results", "Results/Fig. 4", "Results", "Supplementary Fig. 3D", 
    "Supplementary Fig. 3E", "Results/Supplementary Fig. 1", 
    "Results/Fig. 3/Supplementary Table 2", "Results/Fig. 3/Supplementary Table 2", 
    "Results/Fig. 3/Supplementary Table 2", "Results/Fig. 3/Supplementary Table 2", 
    "Results/Fig. 4/Supplementary Table 2", "Results/Fig. 4/Supplementary Table 2", 
    "Results/Fig. 4", "Results/Fig. 4", "Supplementary Table 4", 
    "Supplementary Table 5", "Supplementary Table 4", "Supplementary Table 5", 
    "Extended Data Fig. 4a", "Extended Data Fig. 4a", "Main Fig. 4b / Extended Data Fig. 4b", 
    "Main Fig. 4b / Extended Data Fig. 4b", "Extended Data Fig. 4c", 
    "Extended Data Fig. 4c", "Extended Data Fig. 5a", "Extended Data Fig. 5b", 
    "Results/Fig. 2", "Results/Fig. 3", "Results", "Results", 
    "Fig. 3F", "Results", "Results", "Results", "Results", "Results/Fig. 4A", 
    "Results/Fig. 4B", "Results/Fig. 4D", "Results/Fig. 4E", 
    "Results", "Results", "Results", "Results", "Abstract/Results/Table 2", 
    "Table 2", "Results/Fig. 2B/Table 2", "Supplementary Table S6", 
    "Supplementary Table S5", "Supplementary Table S6", "Fig. 1B / Supplementary Table S7", 
    "Fig. 3", "Figure S10 / Table S10", "Results/Supplementary Fig. S8E", 
    "Results/Supplementary Fig. S8F", "Results/Data Supplement Fig. S2A", 
    "Results/Fig. 3B", "Data Supplement Fig. S2B", "Results/Data Supplement Fig. S3B", 
    "Fig. 4A", "Data Supplement Fig. S4A/S4B", "Table 2", "Table 2", 
    "Table 2", "Table 2", "Table 2", "Table 2", "Table 2", "Table 2", 
    "Table 2", "Supplementary Fig. S1A", "Supplementary Fig. S1B", 
    "Supplementary Fig. S1C", "Supplementary Fig. S1C", "Supplementary Fig. S1C", 
    "Supplementary Fig. S1C", "Supplementary Fig. S2A", "Supplementary Fig. S2A", 
    "Supplementary Fig. S2B", "Supplementary Fig. S2B", "Abstract/Results", 
    "Results", "Abstract/Results", "Results", "Abstract/Results"
    ), Notes = c("Inverted to positive vs negative", "Inverted", 
    "Direction opposite adverse exposure; inverted", "Inverted", 
    "13 relapses; 10 detected before relapse; 2 lacked a sample in the preceding 6 months; 1 detected only after relapse; no adjusted model", 
    "88 positive (18 recurrences); 26 negative (1 recurrence)", 
    "29 positive (7 recurrences); 12 negative (0 recurrences)", 
    "38 positive (6 recurrences); 44 negative (8 recurrences)", 
    "17 positive (6 recurrences); 57 negative (8 recurrences)", 
    "Combined ER+ and TNBC: 9 positive and 36 negative", "All postoperative/follow-up positives relapsed; PPV 100%; median lead time 374 days", 
    "5 positive (3 recurrences); 21 negative (5 recurrences)", 
    "4 positive (3 recurrences); 15 negative (4 recurrences)", 
    "Without threshold: 96 positive (10 distant metastases), 26 negative (0 distant metastases); DRFS HR not estimable", 
    "26 positive (7 distant metastases); 87 negative (5 distant metastases)", 
    "26 positive (7 distant metastases); 87 negative (5 distant metastases)", 
    "20 positive (6 distant metastases); 98 negative (6 distant metastases)", 
    "20 positive (6 distant metastases); 98 negative (6 distant metastases)", 
    "Thresholded high group 53 (9 distant metastases); low group 69 (1 distant metastasis); not binary detection", 
    "Thresholded high group 53 (9 distant metastases); low group 69 (1 distant metastasis)", 
    "Same I-SPY2 source population as Magbanua", "High STB 63; low STB 56", 
    "High STB 63 (11 distant metastases); low STB 56 (1 distant metastasis)", 
    "38 MRD-positive (11 distant metastases); 44 negative (0 distant metastases); DRFS HR not estimable", 
    "Adjusted for pathologic response; non-pCR HR 6.33 (0.82<U+2013>49.12)", 
    "Adjusted for tumor size, nodal status and pCR; MRD remained significant (p=0.008)", 
    "87 cleared; 34 not cleared; standardized as no clearance vs clearance", 
    "89 cleared; 34 not cleared; inverted", "26 positive; 42 negative", 
    "26 positive; 42 negative", "11 always positive; 119 negative at <U+2265>1 timepoint", 
    "11 always positive; 119 negative at <U+2265>1 timepoint", 
    "4 positive; 116 negative; extremely sparse positive group", 
    "4 positive; 116 negative; extremely sparse positive group", 
    "Very sparse; separate from SCANDARE", NA, NA, NA, "Adjusted for pCR, stage, grade, adjuvant therapy", 
    "Adjusted for stage and grade", NA, NA, NA, "28/131 positive", 
    NA, NA, NA, NA, NA, NA, NA, "34/83 positive", "Adjusted for stage, subtype, RCB, age", 
    "Point estimate nearly equals upper CI; likely extraction/source error", 
    "16/83 positive", NA, NA, "44 positive (28 relapses); 20 negative (1 relapse); 29 total relapses", 
    "Inverted", "Adjusted model: RCB3 HR 11.2 (1.32<U+2013>95.64), stage III HR 2.2 (0.83<U+2013>5.82), RCB2 HR 4.6 (0.59<U+2013>35.47), capecitabine HR 0.4 (0.18<U+2013>0.85); global PH test p=0.59", 
    "39 positive (14 relapses); 27 negative (3 relapses); inverted", 
    "39 positive; 27 negative; inverted", "29 T1-positive; 19 events; adjusted for nodal status, pCR, grade and age; AIC 154.31; C-index 0.79", 
    NA, "Adjusted for pCR", NA, "9 positive; 66 negative", "5 positive (2 T1-/T2*+ and 3 T1+/T2*+); 71 negative", 
    "49 positive; 60 negative", NA, NA, "25 positive; 70 negative", 
    NA, NA, "18 positive; 71 negative", NA, NA, "Adjusted for T stage, N stage, grade, MammaPrint, receptor subtype and RCB as applicable", 
    "Adjusted clinicopathologic model", "Full-cohort adjusted model; not the RCB-II/III N=249 analysis", 
    "Full-cohort adjusted model; not the RCB-II/III N=249 analysis", 
    "Wald p=0.213; full-cohort adjusted model", "Wald p=0.269; full-cohort adjusted model", 
    "159 ctDNA-positive patients; p=0.02; total model N not supplied", 
    "81 ctDNA-positive patients; p<0.001; total model N not supplied", 
    "21 ctDNA-positive patients; p=0.001; total model N not supplied", 
    "30 ctDNA-positive patients; p<0.001; total model N not supplied", 
    NA, NA, "Adjusted for pCR", NA, "Adjusted for pCR")), row.names = c(NA, 
-92L), class = "data.frame")
)

# Nonpooled supplementary results retained for the corresponding source-data
# export; these rows are not used in the meta-analytic models.
supplementary_results <- tibble::as_tibble(
structure(list(Study_ID = c("S13", "S13", "S13", "S13", "S13", 
"S13", "S13", "S13", "S13", "S13", "S13", "S13"), Citation = c("Magbanua 2025", 
"Magbanua 2025", "Magbanua 2025", "Magbanua 2025", "Magbanua 2025", 
"Magbanua 2025", "Magbanua 2025", "Magbanua 2025", "Magbanua 2025", 
"Magbanua 2025", "Magbanua 2025", "Magbanua 2025"), `Analysis family` = c("RCB-II/III clearance dynamics", 
"RCB-II/III clearance dynamics", "RCB-II/III clearance dynamics", 
"RCB-II/III clearance dynamics", "RCB-II/III clearance dynamics", 
"Early clearance and RCB-II/III", "Early clearance and RCB-II/III", 
"Early clearance and RCB-II/III", "Early clearance and RCB-II/III", 
"Early clearance and RCB-II/III", "Variant conservation", "Variant conservation"
), `Population/comparison` = c("Persistent ctDNA-negative", "Cleared at T1 (week 3)", 
"Cleared at T2 (week 12)", "Cleared at T3 (post-NAT)", "No clearance", 
"HR+/HER2<U+2212> + paclitaxel", "HR+/HER2<U+2212> + paclitaxel + ICI", 
"TNBC + paclitaxel", "TNBC + paclitaxel + ICI", "HER2+ + anti-HER2 therapy", 
"PT0 vs PT1<U+2032> (N=42)", "PT0 vs PT3 (N=41)"), Outcome = c("DRFS", 
"DRFS", "DRFS", "DRFS", "DRFS", "RCB-II/III outcome", "RCB-II/III outcome", 
"RCB-II/III outcome", "RCB-II/III outcome", "RCB-II/III outcome", 
"Tumor variant conservation", "Tumor variant conservation"), 
    `Effect type` = c("Reference", "HR", "HR", "HR", "HR", "OR", 
    "OR", "OR", "OR", "OR", "Median range", "Median range"), 
    Estimate = c("1", "3.26", "3.18", "11.21", "30.67", "0.39", 
    "0.2", "0.14", "0.29", "0.27", "21.5%<U+2013>28.6%", "18.1%<U+2013>23.3%"
    ), `95% CI / comparator range` = c("NR", "NR", "NR", "NR", 
    "NR", "0.11<U+2013>1.21", "0.04<U+2013>0.82", "0.04<U+2013>0.39", 
    "0.09<U+2013>0.94", "0.05<U+2013>1.26", "90.6%<U+2013>100%", 
    "81.3%<U+2013>100%"), `P value` = c("-", "0.136", "0.151", 
    "<0.001", "<0.001", "0.202", "0.047", "0.001", "0.054", "0.160", 
    "-", "-"), `3-year DRFS (%)` = c(98, 91, 92, 83, 35, NA, 
    NA, NA, NA, NA, NA, NA), `Meta use` = c("No<U+2014>reference", 
    "No<U+2014>CI unavailable", "No<U+2014>CI unavailable", "No<U+2014>CI unavailable", 
    "No<U+2014>CI unavailable", "No<U+2014>pathologic outcome", 
    "No<U+2014>pathologic outcome", "No<U+2014>pathologic outcome", 
    "No<U+2014>pathologic outcome", "No<U+2014>pathologic outcome", 
    "No<U+2014>mechanistic", "No<U+2014>mechanistic"), `Source locator` = c("Supplementary Fig. S3", 
    "Supplementary Fig. S3", "Supplementary Fig. S3", "Supplementary Fig. S3", 
    "Supplementary Fig. S3", "Supplementary Fig. S8", "Supplementary Fig. S8", 
    "Supplementary Fig. S8", "Supplementary Fig. S8", "Supplementary Fig. S8", 
    "Supplementary Table S4/Fig. 6", "Supplementary Table S4/Fig. 6"
    ), `Notes/QC` = c("Complete four-timepoint RCB-II/III analysis N=249", 
    "Reference=persistently negative", "Reference=persistently negative", 
    "Reference=persistently negative", "Reference=persistently negative", 
    "Early clearance at week 3 vs not cleared at week 12", "Not a time-to-event endpoint", 
    "Not a time-to-event endpoint", "Reported CI excludes 1 but p=0.054; retain source-reported values and flag inconsistency", 
    "Not a time-to-event endpoint", "Estimate=all variants; CI field=personalized ctDNA assay variants", 
    "More than 96% of the 16 selected targets were retained overall"
    )), row.names = c(NA, -12L), class = "data.frame")
)
effects <- raw_effects %>%
  transmute(
    effect_id = Effect_ID,
    publication_id = Study_ID,
    study = Citation,
    cohort = Cohort,
    source_timepoint = Timepoint,
    analysis_n = suppressWarnings(as.numeric(`Analysis N`)),
    subgroup_source = Subgroup,
    endpoint = toupper(trimws(as.character(Endpoint))),
    analysis = Analysis,
    reported_contrast = `Reported contrast`,
    hr_reported = as.numeric(HR_reported),
    ci_low_reported = as.numeric(CI_low_reported),
    ci_high_reported = as.numeric(CI_high_reported),
    invert = tolower(trimws(as.character(`Invert?`))),
    cohort_id = `Dependency cluster`,
    dependency_cluster = `Dependency cluster`,
    source_locator = `Source locator`,
    notes = Notes
  ) %>%
  mutate(
    invert = coalesce(invert, "no"),
    hr = if_else(invert == "yes", 1 / hr_reported, hr_reported),
    ci_low = if_else(invert == "yes", 1 / ci_high_reported, ci_low_reported),
    ci_high = if_else(invert == "yes", 1 / ci_low_reported, ci_high_reported),
    yi = log(hr),
    sei = (log(ci_high) - log(ci_low)) / (2 * 1.96),
    vi = sei^2,
    analysis_class = case_when(
      grepl("unadjusted|univariable", analysis, ignore.case = TRUE) ~ "Unadjusted",
      grepl("multivariable|adjusted", analysis, ignore.case = TRUE) ~ "Adjusted",
      TRUE ~ "Unadjusted"
    )
  )
study_dictionary <- tribble(
  ~publication_id, ~year, ~assay_family, ~cohort_subtype,
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
effects <- effects %>% left_join(study_dictionary, by = "publication_id")
# cohort_id is the analytic clustering unit (the underlying patient cohort);
# publication_id/study identify the publication. One cohort may be represented by one or more publications.
cohort_dependency_audit <- effects %>%
  distinct(cohort_id, dependency_cluster, publication_id, study, cohort) %>%
  arrange(cohort_id, publication_id)
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
  "E015", "Postoperative sampling window", "E016", "Postoperative sampling window",
  "E040", "Postoperative sampling window", "E041", "Postoperative sampling window",
  "E047", "Postoperative sampling window", "E048", "Postoperative sampling window",
  "E049", "Postoperative sampling window", "E053", "Postoperative sampling window",
  "E054", "Postoperative sampling window", "E059", "Postoperative sampling window",
  "E060", "Postoperative sampling window", "E081", "Postoperative sampling window",
  "E082", "Postoperative sampling window",
  "E005", "Longitudinal surveillance", "E011", "Longitudinal surveillance",
  "E027", "Longitudinal surveillance", "E042", "Longitudinal surveillance",
  "E043", "Longitudinal surveillance", "E044", "Longitudinal surveillance",
  "E077", "Longitudinal surveillance", "E078", "Longitudinal surveillance",
  "E085", "Longitudinal surveillance", "E086", "Longitudinal surveillance"
)
timepoint_levels <- c(
  "Baseline", "During NAT", "Post-NAT/preoperative",
  "Postoperative sampling window", "Longitudinal surveillance"
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
  # Postoperative sampling window (broad study-defined postoperative timing)
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
  "E041", "E049", "E054", "E082", # postoperative sampling window
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
# Additional analysis restricted to recurrence-specific endpoints with the
# closest clinical definitions across cohorts. Broader DFS, EFS, and IDFS
# composites are excluded because they may include nonrecurrence events.
comparable_recurrence_endpoints <- c("RFS", "RFI", "DRFS", "DRFI", "BCFI")
comparable_recurrence <- primary_timepoint %>%
  filter(toupper(endpoint) %in% comparable_recurrence_endpoints)
if (nrow(comparable_recurrence) != 21L ||
    n_distinct(comparable_recurrence$cohort_id) != 10L) {
  stop("Comparable recurrence set must contain 21 effects from 10 cohorts")
}
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

# Reproducible membership map for every analysis set. Effect selection is fixed
# by Effect_ID and is not performed according to HR magnitude or significance.
tag_analysis_set <- function(d, analysis_set, purpose) {
  d %>%
    transmute(
      analysis_set = analysis_set,
      purpose = purpose,
      effect_id, publication_id, study, cohort_id, source_timepoint,
      clinical_timepoint = as.character(clinical_timepoint), endpoint,
      analysis_class, hr, ci_low, ci_high
    )
}
analysis_set_membership <- bind_rows(
  tag_analysis_set(primary_timepoint, "Primary recurrence-related by sampling window",
                   "One prespecified effect per cohort within each sampling window"),
  tag_analysis_set(comparable_recurrence, "Comparable recurrence-specific endpoints",
                   "Restricted to RFS, RFI, DRFS, DRFI, and BCFI"),
  tag_analysis_set(adjusted_all, "Adjusted estimates by sampling window",
                   "Adjusted estimates pooled only within clinically aligned windows"),
  tag_analysis_set(global_unadjusted, "Independent unadjusted dataset",
                   "One unadjusted estimate per cohort for conventional diagnostics"),
  tag_analysis_set(global_adjusted, "Independent adjusted dataset",
                   "One adjusted estimate per cohort for conventional diagnostics"),
  tag_analysis_set(global_all_studies, "Assay-platform exploratory dataset",
                   "One recurrence-related estimate per underlying cohort"),
  tag_analysis_set(os_primary, "Overall survival exploratory dataset",
                   "One OS estimate per cohort within each sampling window"),
  tag_analysis_set(subtype_data, "Molecular-subtype exploratory dataset",
                   "Subtype estimates analyzed only within clinically aligned windows")
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
  list(primary_timepoint, comparable_recurrence, adjusted_all, global_unadjusted, global_adjusted,
       global_all_studies, os_primary, subtype_data, magbanua_clearance, magbanua_rcb),
  c("primary_timepoint", "comparable_recurrence", "adjusted_all", "global_unadjusted", "global_adjusted",
    "global_all_studies", "os_primary", "subtype_data", "magbanua_clearance", "magbanua_rcb"),
  qc_effects
)
if (anyDuplicated(global_unadjusted$cohort_id)) stop("global_unadjusted is not independent by cohort")
if (anyDuplicated(global_adjusted$cohort_id)) stop("global_adjusted is not independent by cohort")
if (anyDuplicated(global_all_studies$cohort_id)) stop("global_all_studies is not independent by cohort")
cohort_cluster_map <- cohort_dependency_audit %>%
  distinct(cohort_id, dependency_cluster)
if (any(is.na(cohort_cluster_map$cohort_id) | is.na(cohort_cluster_map$dependency_cluster))) {
  stop("Missing cohort_id or dependency_cluster in cohort audit")
}
if (any(cohort_cluster_map$cohort_id != cohort_cluster_map$dependency_cluster)) {
  stop("cohort_id and dependency_cluster disagree")
}
dup_primary <- primary_timepoint %>% count(cohort_id, clinical_timepoint) %>% filter(n > 1)
if (nrow(dup_primary)) stop("Primary timepoint set has duplicate cohort-timepoint rows")
# Endpoint-selection audit for Supplementary eTable 4. Candidate outcomes are
# listed within each underlying cohort and clinical sampling window.
candidate_outcomes <- effects %>%
  filter(!is.na(clinical_timepoint), !is.na(endpoint)) %>%
  group_by(cohort_id, clinical_timepoint) %>%
  summarise(
    reported_candidate_outcomes = paste(sort(unique(endpoint)), collapse = ", "),
    .groups = "drop"
  )
cohort_audit_counts <- effects %>%
  group_by(cohort_id) %>%
  summarise(
    publications_in_cohort = n_distinct(publication_id),
    .groups = "drop"
  )
selected_window_counts <- primary_timepoint %>%
  count(cohort_id, name = "selected_effects_from_cohort")
endpoint_selection_audit <- primary_timepoint %>%
  left_join(candidate_outcomes, by = c("cohort_id", "clinical_timepoint")) %>%
  left_join(cohort_audit_counts %>% select(cohort_id, publications_in_cohort), by = "cohort_id") %>%
  left_join(selected_window_counts, by = "cohort_id") %>%
  mutate(
    harmonized_comparison = "Unfavorable ctDNA state vs favorable ctDNA state",
    selection_rationale = case_when(
      analysis_class == "Adjusted" ~ "Adjusted estimate preferred for the same cohort-window contrast.",
      TRUE ~ "Direct binary ctDNA contrast; study-designated or most complete recurrence-related endpoint."
    ),
    comparable_recurrence_status = if_else(
      toupper(endpoint) %in% comparable_recurrence_endpoints,
      "Included in comparable recurrence-specific analysis",
      "Excluded from comparable analysis because the endpoint was a broader composite"
    ),
    same_underlying_cohort = case_when(
      publications_in_cohort > 1 ~ "Yes; cohort represented in multiple publications",
      selected_effects_from_cohort > 1 ~ "Yes; serial effects from the same cohort",
      TRUE ~ "No other selected effect from this cohort"
    )
  ) %>%
  transmute(
    publication_id, publication = study, underlying_cohort_id = cohort_id,
    cohort_description = cohort, effect_id, sampling_window = clinical_timepoint,
    reported_candidate_outcomes, selected_endpoint = endpoint,
    selected_hr = hr, selected_ci_low = ci_low, selected_ci_high = ci_high,
    adjustment = analysis_class, reported_comparison = reported_contrast,
    harmonized_comparison, selection_rationale, comparable_recurrence_status,
    same_underlying_cohort
  )
write.csv(effects, file.path(out_dir, "Source_Data_all_standardized_effects.csv"), row.names = FALSE)
write.csv(cohort_dependency_audit,
          file.path(out_dir, "Cohort_dependency_cluster_audit.csv"), row.names = FALSE)
write.csv(endpoint_selection_audit,
          file.path(out_dir, "eTable_2_endpoint_HR_selection_audit.csv"), row.names = FALSE)
write.csv(analysis_set_membership,
          file.path(out_dir, "Analysis_set_membership_audit.csv"), row.names = FALSE)
write.csv(primary_timepoint, file.path(out_dir, "Source_Data_primary_timepoint.csv"), row.names = FALSE)
write.csv(comparable_recurrence,
          file.path(out_dir, "Source_Data_comparable_recurrence.csv"), row.names = FALSE)
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
  for (s in unique(d$cohort_id)) {
    idx <- which(d$cohort_id == s)
    if (length(idx) > 1) {
      V[idx, idx] <- rho * sqrt(outer(d$vi[idx], d$vi[idx]))
      diag(V)[idx] <- d$vi[idx]
    }
  }
  V
}
fit_continuum_rve <- function(d, rho = 0.60, moderator = TRUE) {
  d <- d %>% arrange(cohort_id, clinical_timepoint, effect_id)
  V <- make_sampling_V(d, rho)
  if (moderator) {
    # The model matrix already contains one coefficient per timepoint, so the
    # default intercept must be disabled to avoid redundant predictors.
    mods <- model.matrix(~ 0 + clinical_timepoint, data = d)
    fit <- metafor::rma.mv(
      yi = d$yi, V = V, mods = mods, intercept = FALSE,
      random = ~ 1 | cohort_id/effect_id,
      method = "REML", data = d
    )
  } else {
    # rma.mv does not accept an explicitly supplied mods=NULL object. Omit the
    # argument entirely for the intercept-only global model.
    fit <- metafor::rma.mv(
      yi = d$yi, V = V,
      random = ~ 1 | cohort_id/effect_id,
      method = "REML", data = d
    )
  }
  robust <- clubSandwich::coef_test(
    fit, vcov = "CR2", cluster = d$cohort_id, test = "Satterthwaite"
  )
  list(data = d, V = V, fit = fit, robust = robust, rho = rho)
}
continuum_rve <- fit_continuum_rve(primary_timepoint, ASSUMED_WITHIN_COHORT_RHO, TRUE)
continuum_global_rve <- fit_continuum_rve(primary_timepoint, ASSUMED_WITHIN_COHORT_RHO, FALSE)
comparable_rve <- fit_continuum_rve(comparable_recurrence, ASSUMED_WITHIN_COHORT_RHO, TRUE)
comparable_global_rve <- fit_continuum_rve(comparable_recurrence, ASSUMED_WITHIN_COHORT_RHO, FALSE)
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
          file.path(out_dir, "Sensitivity_within_cohort_rho.csv"), row.names = FALSE)
comparable_timepoint_results <- fit_by_group(
  comparable_recurrence, clinical_timepoint,
  "Comparable recurrence-specific endpoints by timepoint"
)
comparable_global_ct <- as.data.frame(comparable_global_rve$robust)
comparable_global_critical_t <- qt(0.975, df = comparable_global_ct$df_Satt[1])
comparable_global_result <- tibble(
  model = "Comparable recurrence-specific endpoints across timepoints",
  group = "Overall",
  k = nrow(comparable_recurrence),
  cohorts = n_distinct(comparable_recurrence$cohort_id),
  log_hr = comparable_global_ct$beta[1],
  hr = exp(comparable_global_ct$beta[1]),
  ci_low = exp(comparable_global_ct$beta[1] - comparable_global_critical_t * comparable_global_ct$SE[1]),
  ci_high = exp(comparable_global_ct$beta[1] + comparable_global_critical_t * comparable_global_ct$SE[1]),
  p_value = comparable_global_ct$p_Satt[1],
  tau2 = NA_real_, i2 = NA_real_, q_p_value = NA_real_,
  method = "Multilevel REML + CR2; cohort clustered; rho=0.60"
)
comparable_etable5 <- bind_rows(
  comparable_timepoint_results %>%
    mutate(effects = k, cohorts = k) %>%
    select(analysis = group, effects, cohorts, hr, ci_low, ci_high, p_value, i2),
  comparable_global_result %>%
    transmute(
      analysis = "Across time points (additional CR2 analysis)",
      effects = k, cohorts, hr, ci_low, ci_high, p_value, i2
    )
)
write.csv(comparable_etable5,
          file.path(out_dir, "eTable_4_comparable_recurrence_results.csv"), row.names = FALSE)
write.csv(as.data.frame(comparable_global_rve$robust),
          file.path(out_dir, "Model_comparable_recurrence_global_CR2.csv"), row.names = FALSE)
fit_uni_global <- fit_reml(global_unadjusted)
fit_adj_global <- fit_reml(global_adjusted)
model_results <- bind_rows(
  fit_by_group(primary_timepoint, clinical_timepoint, "Primary recurrence-related by timepoint"),
  comparable_timepoint_results,
  comparable_global_result,
  fit_by_group(os_primary, clinical_timepoint, "OS by timepoint"),
  fit_by_group(subtype_data, subtype_label, "Molecular subtype exploratory"),
  fit_by_group(global_all_studies, assay_family, "Assay platform exploratory"),
  model_row(fit_uni_global, "Independent unadjusted diagnostic set"),
  model_row(fit_adj_global, "Independent adjusted diagnostic set")
)

# Locked-result checks corresponding to the submitted manuscript. These checks
# stop execution if the workbook or effect-selection vectors change silently.
assert_close <- function(actual, expected, label, tolerance = 0.01) {
  if (length(actual) != 1L || !is.finite(actual) || abs(actual - expected) > tolerance) {
    stop(label, " does not match the locked manuscript value. Expected ",
         expected, "; obtained ", paste(actual, collapse = ", "))
  }
}
get_hr <- function(model_name, group_name) {
  model_results %>%
    filter(model == model_name, group == group_name) %>%
    pull(hr)
}
assert_close(
  get_hr("Primary recurrence-related by timepoint", "Postoperative sampling window"),
  6.45, "Postoperative recurrence-related HR"
)
assert_close(
  get_hr("Primary recurrence-related by timepoint", "Longitudinal surveillance"),
  28.31, "Surveillance recurrence-related HR"
)
assert_close(
  get_hr("OS by timepoint", "Postoperative sampling window"),
  11.30, "Postoperative OS HR"
)
assert_close(rho_sensitivity$hr[rho_sensitivity$rho == 0.60],
             6.65, "Cross-window CR2 HR")
assert_close(comparable_global_result$hr, 5.91,
             "Comparable-endpoint cross-window HR")
assert_close(exp(as.numeric(fit_uni_global$b)), 13.30,
             "Independent unadjusted HR")
assert_close(exp(as.numeric(fit_adj_global$b)), 5.98,
             "Independent adjusted HR")
write.csv(model_results, file.path(out_dir, "Meta_model_results.csv"), row.names = FALSE)
# ------------------------------ 5. Forest-plot helpers ------------------------
format_p_value <- function(p) {
  if (length(p) == 0 || is.na(p) || !is.finite(p)) return("P=NA")
  if (p < 0.001) return("P<0.001")
  sprintf("P=%.3f", p)
}
format_bound <- function(x) {
  ifelse(
    x < 0.01 | x >= 1000,
    format(x, digits = 3, scientific = TRUE, trim = TRUE),
    sprintf("%.2f", x)
  )
}
format_hr_ci_p <- function(hr, ci_low, ci_high, p) {
  sprintf("%.2f (%s-%s); %s", hr, format_bound(ci_low), format_bound(ci_high),
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
          publication_id = NA_character_, cohort_id = NA_character_, study = NA_character_, endpoint = NA_character_,
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
  "Postoperative sampling window" = unname(palette_ctdna["landmark"]),
  "Longitudinal surveillance" = unname(palette_ctdna["surveillance"])
)
# ------------------------------ 6. Main figures -------------------------------
fig1a <- forest_plot(
  primary_timepoint, clinical_timepoint,
  "Recurrence-related outcomes across clinical sampling windows",
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
# The subtype analysis is supplementary; panel-specific exports are created below.
fig3b <- forest_plot(
  global_all_studies, assay_family,
  "Assay-platform subgroup analysis",
  "Exploratory and potentially confounded by timepoint, cohort and analysis type",
  pool = TRUE,
  colour_map = c("Sequencing-based" = "#3B6FB6", "dPCR-based" = "#D98E3D")
)
save_pub(fig3b, "Figure_3_Assay_Forest", 230, 125)
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
fig4b <- funnel_plot(
  global_adjusted, fit_adj_global, "Funnel plot: independent adjusted estimates",
  if (is.null(egger_adj)) NA_real_ else egger_adj$p_value
)
save_pub_pair(fig4a, fig4b, "Supplementary_Figure_15_Funnel_Plots", 230, 105)
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
  save_pub(p_egger, "Supplementary_Figure_16_Egger_Regression", 135, 105)
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
           "Studies added chronologically; REML + Hartung-Knapp after the first study.\n",
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
         "Supplementary_Figure_12_Cumulative_Unadjusted", 215, 120)
save_pub(plot_cumulative(cum_adj, "Cumulative meta-analysis: adjusted estimates"),
         "Supplementary_Figure_11_Cumulative_Adjusted", 210, 95)
# ------------------------------ 9. Influence diagnostics ----------------------
save_base_plot(
  function() {
    par(mar = c(5.0, 4.2, 2.0, 1.0), family = "sans")
    metafor::baujat(fit_uni_global, symbol = "ids", cex = 0.85,
                    xlab = expression(Delta*Q), ylab = expression(Influence~on~pooled~effect))
    title("Baujat plot: unadjusted independent set",
          sub = fit_annotation(fit_uni_global), cex.main = 0.95, cex.sub = 0.70)
  },
  "Supplementary_Figure_13A_Baujat_Unadjusted", 120, 105
)
save_base_plot(
  function() {
    par(mar = c(5.0, 4.2, 2.0, 1.0), family = "sans")
    metafor::baujat(fit_adj_global, symbol = "ids", cex = 0.85,
                    xlab = expression(Delta*Q), ylab = expression(Influence~on~pooled~effect))
    title("Baujat plot: adjusted independent set",
          sub = fit_annotation(fit_adj_global), cex.main = 0.95, cex.sub = 0.70)
  },
  "Supplementary_Figure_13B_Baujat_Adjusted", 120, 105
)
save_base_plot(
  function() {
    par(mar = c(5.0, 4.2, 2.0, 1.0), family = "sans")
    metafor::radial(fit_uni_global)
    title("Radial plot: unadjusted independent set",
          sub = fit_annotation(fit_uni_global), cex.main = 0.95, cex.sub = 0.70)
  },
  "Supplementary_Figure_5A_Radial_Unadjusted", 120, 105
)
save_base_plot(
  function() {
    par(mar = c(5.0, 4.2, 2.0, 1.0), family = "sans")
    metafor::radial(fit_adj_global)
    title("Radial plot: adjusted independent set",
          sub = fit_annotation(fit_adj_global), cex.main = 0.95, cex.sub = 0.70)
  },
  "Supplementary_Figure_5B_Radial_Adjusted", 120, 105
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
save_pub(p_loo, "Supplementary_Figure_10_Leave_One_Out", 210, 115)
# ------------------------------ 10. Timepoint-specific supplementary forests --
timepoint_file_tags <- c(
  "Baseline" = "3A_Baseline",
  "During NAT" = "3B_During_NAT",
  "Post-NAT/preoperative" = "3C_Preoperative",
  "Postoperative sampling window" = "3D_Postoperative_Sampling_Window",
  "Longitudinal surveillance" = "3E_Surveillance"
)
for (tp in timepoint_levels) {
  d_tp <- primary_timepoint %>% filter(clinical_timepoint == tp) %>% mutate(panel = tp)
  if (nrow(d_tp)) {
    p_tp <- forest_plot(
      d_tp, panel, paste0(tp, " ctDNA"),
      "Recurrence-related time-to-event outcomes", pool = TRUE,
      colour_map = setNames(timepoint_colours[tp], tp)
    )
    save_pub(p_tp, paste0("Supplementary_Figure_", timepoint_file_tags[tp]), 220,
             max(70, 45 + 8 * nrow(d_tp)))
  }
}
# ------------------------------ 11. Dependency-corrected subtype panels ------
# HER2+ and HR+/HER2- estimates are displayed descriptively without pooling.
# TNBC estimates are pooled only within sampling windows, where each underlying
# patient cohort contributes no more than one estimate. No across-time-point
# conventional subtype summary is calculated.
ds <- subtype_data
p_jama <- function(p) {
  if (is.na(p)) return("P = NA")
  if (p < .001) return("P < .001")
  paste0("P = ", sub("^0", "", sprintf("%.3f", p)))
}
format_upper <- format_bound
append_timepoint_pools <- function(dat) {
  preferred_order <- c("During NAT", "Post-NAT/preoperative", "Postoperative sampling window")
  dat$clinical_timepoint <- factor(dat$clinical_timepoint,
                                   levels = preferred_order)
  dat <- dat %>% arrange(clinical_timepoint)
  out <- lapply(preferred_order, function(tp) {
    g <- dat %>% filter(clinical_timepoint == tp)
    if (!nrow(g)) return(NULL)
    # Ordinary pooling is valid within a sampling window only when every
    # underlying cohort contributes no more than one estimate.
    cluster_var <- if ("cohort_id" %in% names(g)) "cohort_id" else "dependency_cluster"
    if (anyDuplicated(g[[cluster_var]])) {
      stop("Duplicate underlying cohort within TNBC sampling window: ", tp)
    }
    fit <- metafor::rma.uni(
      yi = g$yi, vi = g$vi,
      method = "REML", test = "knha"
    )
    pooled <- g[1, , drop = FALSE]
    pooled[,] <- NA
    pooled$clinical_timepoint <- tp
    pooled$hr <- exp(as.numeric(fit$b))
    pooled$ci_low <- exp(fit$ci.lb)
    pooled$ci_high <- exp(fit$ci.ub)
    pooled$yi <- as.numeric(fit$b)
    pooled$sei <- as.numeric(fit$se)
    pooled$vi <- as.numeric(fit$se)^2
    pooled$p_value_panel <- fit$pval
    pooled$pooled <- TRUE
    pooled$display <- paste0(
      tp, " random effects (k=", fit$k,
      "; I\u00B2=", sprintf("%.1f", fit$I2), "%)"
    )
    g$pooled <- FALSE
    g$p_value_panel <- 2 * pnorm(-abs(g$yi / g$sei))
    g$display <- paste0(g$study, " (", g$endpoint, "; ", g$analysis_class, ")")
    bind_rows(g, pooled)
  })
  bind_rows(out)
}
prepare_descriptive <- function(dat) {
  dat %>%
    mutate(
      pooled = FALSE,
      p_value_panel = 2 * pnorm(-abs(yi / sei)),
      display = paste0(study, " (", endpoint, "; ", analysis_class, ")")
    )
}
forest_panel_corrected <- function(dat, label, tag, limits, breaks,
                                   subtitle = NULL) {
  dat$stat <- paste0(
    sprintf("%.2f (%s-%s); ", dat$hr, format_bound(dat$ci_low),
            format_bound(dat$ci_high)),
    vapply(dat$p_value_panel, p_jama, character(1))
  )
  dat$row <- factor(seq_len(nrow(dat)), levels = rev(seq_len(nrow(dat))))
  pm <- ggplot(dat, aes(hr, row)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "#777777") +
    geom_segment(
      aes(x = pmax(ci_low, limits[1]),
          xend = pmin(ci_high, limits[2]), yend = row),
      colour = "#3B6FB6", linewidth = .6
    ) +
    geom_point(aes(shape = pooled), colour = "#3B6FB6", size = 3) +
    scale_shape_manual(values = c(`FALSE` = 15, `TRUE` = 18), guide = "none") +
    scale_x_log10(limits = limits, breaks = breaks,
                  labels = label_number()) +
    scale_y_discrete(labels = setNames(dat$display, seq_len(nrow(dat)))) +
    labs(title = label, subtitle = subtitle, tag = tag,
         x = "Hazard ratio (95% CI)", y = NULL) +
    theme_classic(base_family = "Arial", base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 10, lineheight = 1.0,
                                   margin = margin(b = 5)),
      plot.tag = element_text(family = "Arial", face = "bold", size = 10),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 10),
      plot.margin = margin(6, 5, 6, 7)
    )
  pt <- ggplot(dat, aes(0, row, label = stat,
                        fontface = ifelse(pooled, "bold", "plain"))) +
    geom_text(hjust = 0, size = 3.5) +
    scale_y_discrete() +
    scale_x_continuous(limits = c(0, 2.5)) +
    theme_void(base_family = "Arial") +
    theme(plot.margin = margin(18, 5, 10, 2))
  pm + pt + plot_layout(widths = c(2.3, 1.35))
}
d_a <- ds %>% filter(subtype_label == "HER2+") %>% prepare_descriptive()
d_b <- ds %>% filter(subtype_label == "HR+/HER2-") %>% prepare_descriptive()
d_c <- ds %>% filter(subtype_label == "TNBC") %>% append_timepoint_pools()
p11a <- forest_panel_corrected(
  d_a, "HER2+ subtype", "A",
  c(.1, 1e4), c(.1, 1, 10, 100, 1000, 10000),
  "Two estimates from a single underlying cohort are shown descriptively and were not pooled."
)
p11b <- forest_panel_corrected(
  d_b, "HR+/HER2- subtype", "B",
  c(.5, 100), c(.5, 1, 5, 25, 100)
)
p11c <- forest_panel_corrected(
  d_c, "TNBC subtype", "C",
  c(.5, 2000), c(.5, 1, 10, 100, 1000)
)
save_pub(p11a, "Supplementary_Figure_9A_HER2", 320, 85)
save_pub(p11b, "Supplementary_Figure_9B_HR_HER2", 320, 65)
save_pub(p11c, "Supplementary_Figure_9C_TNBC", 360, 145)
# Magbanua 2025: correlated estimates from one cohort are displayed but not pooled.
p_mag_clearance <- forest_plot(
  magbanua_clearance, panel,
  "Magbanua 2025: adjusted ctDNA clearance trajectory",
  "Common reference: persistently ctDNA-negative; correlated estimates are not meta-analysed",
  pool = FALSE,
  colour_map = c("Full-cohort adjusted clearance trajectory" = "#8A5AA5")
)
save_pub(p_mag_clearance, "Supplementary_Figure_7_Magbanua_Clearance_Trajectory", 220, 85)
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
save_pub(p_mag_rcb, "Supplementary_Figure_8_Magbanua_RCB_Strata", 220, 90)
# ------------------------------ 11. Trim-and-fill -----------------------------
tf_uni <- metafor::trimfill(fit_uni_global)
small_study_effects <- tibble(
  analysis = c(
    "Independent recurrence-related, unadjusted",
    "Independent recurrence-related, adjusted",
    "Overall survival"
  ),
  k = c(nrow(global_unadjusted), nrow(global_adjusted), n_distinct(os_primary$cohort_id)),
  egger_intercept = c(
    if (is.null(egger_uni)) NA_real_ else egger_uni$intercept,
    if (is.null(egger_adj)) NA_real_ else egger_adj$intercept,
    NA_real_
  ),
  egger_p = c(
    if (is.null(egger_uni)) NA_real_ else egger_uni$p_value,
    if (is.null(egger_adj)) NA_real_ else egger_adj$p_value,
    NA_real_
  ),
  formal_test = c(
    if (is.null(egger_uni)) "Not performed (k < 10)" else "Performed",
    if (is.null(egger_adj)) "Not performed (k < 10)" else "Performed",
    "Not performed (k < 10)"
  ),
  trimfill_missing = c(tf_uni$k0, NA_integer_, NA_integer_),
  original_hr = c(exp(as.numeric(fit_uni_global$b)), exp(as.numeric(fit_adj_global$b)), NA_real_),
  original_ci_low = c(exp(fit_uni_global$ci.lb), exp(fit_adj_global$ci.lb), NA_real_),
  original_ci_high = c(exp(fit_uni_global$ci.ub), exp(fit_adj_global$ci.ub), NA_real_),
  trimfill_hr = c(exp(as.numeric(tf_uni$b)), NA_real_, NA_real_),
  trimfill_ci_low = c(exp(tf_uni$ci.lb), NA_real_, NA_real_),
  trimfill_ci_high = c(exp(tf_uni$ci.ub), NA_real_, NA_real_)
)
write.csv(small_study_effects,
          file.path(out_dir, "eTable_5_small_study_effects.csv"), row.names = FALSE)
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
  "Supplementary_Figure_17_Trim_and_Fill", 230, 95
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
save_pub(p_weights, "Supplementary_Figure_6_Study_Weights", 155, 110)
rve_display <- rho_sensitivity %>%
  slice(which.min(abs(rho - ASSUMED_WITHIN_COHORT_RHO)))
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
    "Postoperative sampling window" = "Postoperative\nsampling window",
    "Longitudinal surveillance" = "Longitudinal\nsurveillance"
  )) +
  scale_colour_manual(values = timepoint_colours, guide = "none") +
  labs(x = NULL, y = "Hazard ratio (log scale)",
       title = "Distribution of effect estimates by ctDNA timepoint",
       subtitle = distribution_subtitle) +
  theme_nature(base_size = 7) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, lineheight = 0.9))
save_pub(p_distribution, "Supplementary_Figure_4_Timepoint_Distribution", 183, 105)
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
  draws_raw <- posterior::as_draws_df(bayes_fit)
  draws <- draws_raw %>%
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
  bayes_diagnostics <- posterior::summarise_draws(
    posterior::subset_draws(draws_raw, variable = c("b_Intercept", "sigma")),
    "mean", "sd", "rhat", "ess_bulk", "ess_tail"
  )
  nuts <- brms::nuts_params(bayes_fit)
  sampler_diagnostics <- tibble(
    divergent_transitions = sum(
      nuts$Parameter == "divergent__" & nuts$Value == 1
    ),
    maximum_treedepth_hits = sum(
      nuts$Parameter == "treedepth__" & nuts$Value >= 12
    ),
    post_warmup_draws = posterior::ndraws(draws_raw)
  )
  write.csv(bayes_diagnostics,
            file.path(out_dir, "Bayesian_convergence_diagnostics.csv"), row.names = FALSE)
  write.csv(sampler_diagnostics,
            file.path(out_dir, "Bayesian_sampler_diagnostics.csv"), row.names = FALSE)
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
  save_pub(p_bayes, "Supplementary_Figure_14_Bayesian_Posterior", 135, 100)
} else {
  if (!isTRUE(RUN_BAYESIAN)) {
    message("Bayesian analysis skipped because CTDNA_META_RUN_BAYESIAN=false.")
  } else {
    message(
      "Bayesian analysis skipped because brms and/or posterior could not be loaded. ",
      "Install their dependencies and rerun with CTDNA_META_RUN_BAYESIAN=true."
    )
  }
}
# ------------------------------ 14. PRISMA flow diagram -----------------------
# Final PRISMA counts reported in the manuscript and Supplement 1.
prisma_counts <- tibble(
  identified_databases = 1183L,
  identified_other = 32L,
  duplicates_removed = 344L,
  records_screened = 871L,
  records_excluded = 768L,
  reports_sought = 103L,
  reports_not_retrieved = 32L,
  reports_assessed = 71L,
  reports_excluded = 57L,
  studies_included = 14L
)
write.csv(prisma_counts, file.path(out_dir, "PRISMA_counts.csv"), row.names = FALSE)
# The PRISMA diagram is a manuscript-layout artifact rather than a statistical
# analysis output. Its verified counts are exported above; the diagram itself is
# maintained in Supplement 1 and is not regenerated by this script.
# ------------------------------ 15. Audit trail and completion ----------------
capture.output(
  list(
    primary_timepoint_models = fit_by_group(
      primary_timepoint, clinical_timepoint, "Primary recurrence-related by timepoint"
    ),
    comparable_recurrence_timepoint_models = comparable_timepoint_results,
    comparable_recurrence_global_CR2 = comparable_global_rve$robust,
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
cat("Embedded effect-level records: ", nrow(raw_effects), "\n", sep = "")
cat("Output directory: ", out_dir, "\n", sep = "")
cat("Important: across-timepoint conventional pooling was replaced by multilevel CR2 analysis.\n")
cat("Important: Egger testing was restricted to an independent set with k >= 10.\n")
cat("Figure export enabled: ", MAKE_FIGURES, "\n", sep = "")

