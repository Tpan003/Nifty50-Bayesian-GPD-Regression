# ==============================================================================
# FORECASTING EXTREME NIFTY 50 LOSSES — COMPLETE SEQUENTIAL R PIPELINE
# ==============================================================================
# This file combines the entire reproducible project in exact execution order:
#   Step 0: Package Verification
#   Step 1: Data Import, Date Validation & Lagged Covariate Construction (R/data_functions.R)
#   Step 2: GPD Likelihood, Priors, MAP Estimation & Rolling Evaluation (R/model_functions.R)
#   Step 3: Shared Project Configuration (R/project_config.R)
#   Step 4: Independent Runtime Validation Suite (R/validation_functions.R)
#   Step 5: Prespecified Reduced-Model Comparison A-D (R/reduced_model_comparison.R)
#   Step 6: Interactive File/Folder Pop-Ups & Main Pipeline (analysis.R)
#   Step 7: Predictor Multicollinearity Audit (R/collinearity_diagnostics.R)
#   Step 8: Presentation Tables & Publication Plots (R/create_presentation_materials.R)
#   Step 9: Report Summary Tables & Diagnostic Figures (report.Rmd chunks)
#   Step 10: Report Tables 2-6, 8, conditional probabilities and Figures 2-15 (NEW)
# Audit fixes applied 2026-11-02 are tagged [FIX B-xx] (see audit log).
# ==============================================================================
