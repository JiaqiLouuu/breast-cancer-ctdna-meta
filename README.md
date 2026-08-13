# breast-cancer-ctdna-meta
Meta analysis of ctDNA in early breast cancer
# Tumor-informed ctDNA Meta-Analysis in Early Breast Cancer

Meta-analysis of tumor-informed circulating tumor DNA (ctDNA) across the neoadjuvant continuum in early breast cancer.

## Overview

This repository contains the reproducible R pipeline for a corrected meta-analysis and publication-figure workflow. The analysis evaluates the prognostic value of unfavourable tumor-informed ctDNA status on time-to-event outcomes, with repeated timepoints handled via multilevel random-effects models.

## Contents

- `main_analysis.R` — Core analysis script covering:
  - Timepoint-specific forest plots (baseline, during NAT, post-NAT/preoperative, post-surgery landmark, longitudinal surveillance)
  - Adjusted and unadjusted hazard-ratio pooling (REML + Hartung-Knapp)
  - Multilevel robust variance estimation (RVE) with CR2 small-sample correction
  - Overall survival, molecular subtype, and assay-platform exploratory analyses
  - Funnel plots, Egger regression, trim-and-fill, and cumulative meta-analysis
  - Influence diagnostics (Baujat, radial, leave-one-out)
  - Bayesian random-effects sensitivity analysis (brms)
  - PRISMA flow diagram template

## Data Requirements

The script expects an Excel workbook named `tumor_informed_ctDNA_meta_extraction_14_studies.xlsx` containing:
- `Effect_Estimates` sheet — extracted HRs with 95% CIs
- `Supplementary_Results` sheet — supplementary data

Place the workbook in the same directory as the script before running.

## Dependencies

- **R** ≥ 4.3
- **Required packages**: `metafor`, `clubSandwich`, `ggplot2`, `ggrepel`, `dplyr`, `tidyr`, `readxl`, `svglite`, `ragg`, `scales`
- **Optional packages**: `brms`, `posterior` (for Bayesian sensitivity analysis), `DiagrammeR` (for PRISMA diagram)

## Key Statistical Methods

- **Primary model**: REML random effects with Hartung-Knapp (Knapp-Hartung) inference
- **Repeated timepoints**: Multilevel model (`rma.mv`) + CR2 robust standard errors (`clubSandwich`)
- **Assumed within-study correlation**: ρ = 0.60 (sensitivity: 0, 0.30, 0.60, 0.90)
- **Heterogeneity**: τ² and I² reported for all pooled estimates

## Outputs

Running the script generates:
- Publication-ready figures (SVG, PDF, TIFF, PNG): forest plots, funnel plots, cumulative analyses, diagnostic plots
- Source data CSVs for all standardized effect sets
- Model audit files and session info for reproducibility

## Citation

If you use this code, please cite the accompanying publication.

## Contact

For questions regarding the analysis, please open an issue or contact the corresponding author.
