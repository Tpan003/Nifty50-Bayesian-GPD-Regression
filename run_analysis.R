# Main entry point for the Nifty 50 extreme-loss project.
# Run from the repository root:
#   Rscript run_analysis.R "data/RT_ONLY_COVARIATES.xlsx" "."

required_packages <- c("readxl", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running the project: ",
    paste(missing_packages, collapse = ", ")
  )
}

# Load the project functions in dependency order.
function_files <- c(
  "R/model_functions.R",
  "R/project_config.R",
  "R/data_functions.R",
  "R/validation_functions.R",
  "R/reduced_model_comparison.R",
  "R/collinearity_diagnostics.R",
  "R/presentation_materials.R",
  "R/report_materials.R"
)

missing_files <- function_files[!file.exists(function_files)]
if (length(missing_files) > 0L) {
  stop(
    "Run this script from the project root. Missing files: ",
    paste(missing_files, collapse = ", ")
  )
}
invisible(lapply(function_files, source, local = .GlobalEnv))

select_input_file <- function() {
  message("Choose the input Excel workbook.")

  if (exists("choose.files", envir = asNamespace("utils"), mode = "function")) {
    path <- utils::choose.files(
      caption = "Select the Nifty 50 input workbook",
      multi = FALSE,
      filters = matrix(
        c("Excel workbooks", "*.xlsx"),
        ncol = 2,
        byrow = TRUE
      )
    )
  } else {
    path <- utils::file.choose()
  }

  if (length(path) == 0L || is.na(path) || !nzchar(path)) {
    stop("No input workbook was selected.")
  }
  path
}

select_output_folder <- function(default_path = getwd()) {
  if (exists("choose.dir", envir = asNamespace("utils"), mode = "function")) {
    path <- utils::choose.dir(
      default = default_path,
      caption = "Choose the folder for analysis outputs"
    )
  } else if (requireNamespace("tcltk", quietly = TRUE)) {
    path <- tcltk::tk_choose.dir(
      default = default_path,
      caption = "Choose the folder for analysis outputs"
    )
  } else {
    message("Folder selection is unavailable. Outputs will be saved in: ", default_path)
    path <- default_path
  }

  if (is.null(path) || length(path) == 0L || is.na(path) || !nzchar(path)) {
    stop("No output folder was selected.")
  }
  path
}

arguments <- commandArgs(trailingOnly = TRUE)
input_from_environment <- Sys.getenv("NIFTY_INPUT", unset = "")
output_from_environment <- Sys.getenv("NIFTY_OUTPUT", unset = "")

input_file <- if (length(arguments) >= 1L) {
  arguments[1]
} else if (nzchar(input_from_environment)) {
  input_from_environment
} else if (interactive()) {
  select_input_file()
} else {
  stop(
    "Provide the workbook path as the first argument or set NIFTY_INPUT.\n",
    "Example: Rscript run_analysis.R data/RT_ONLY_COVARIATES.xlsx ."
  )
}

input_file <- normalizePath(input_file, winslash = "/", mustWork = TRUE)

output_directory <- if (length(arguments) >= 2L) {
  arguments[2]
} else if (nzchar(output_from_environment)) {
  output_from_environment
} else if (interactive()) {
  select_output_folder(dirname(input_file))
} else {
  getwd()
}

dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)
output_directory <- normalizePath(
  output_directory, winslash = "/", mustWork = TRUE
)

results_directory <- file.path(output_directory, "results")
dir.create(results_directory, showWarnings = FALSE, recursive = TRUE)

write_project_csv <- function(data, path) {
  data_to_write <- data
  if ("Date" %in% names(data_to_write)) {
    data_to_write$Date <- format(as.Date(data_to_write$Date), "%Y-%m-%d")
  }
  utils::write.csv(
    data_to_write,
    path,
    row.names = FALSE,
    na = ""
  )
}

message("Input workbook: ", input_file)
message("Output folder:  ", output_directory)

# Step 1: prepare the data and create strictly lagged predictor columns.
config <- project_config()
merged_data <- prepare_and_merge(
  input_file,
  threshold = config$threshold
)

# Step 2: estimate the baseline and prior-regularised GPD models on rolling
# training windows, then score the held-out exceedance observations.
set.seed(config$seed)
predictions <- rolling_evaluation(merged_data, config)
metrics <- score_predictions(predictions)
fit_diagnostics <- attr(predictions, "fit_diagnostics")

# Save the core results before validation so that a failed check can be traced.
write_project_csv(
  merged_data,
  file.path(results_directory, "merged_data.csv")
)
write_project_csv(
  build_extremes(merged_data),
  file.path(results_directory, "extreme_severity_data.csv")
)
write_project_csv(
  predictions,
  file.path(results_directory, "rolling_predictions.csv")
)
write_project_csv(
  metrics,
  file.path(results_directory, "rolling_metrics.csv")
)
write_project_csv(
  fit_diagnostics,
  file.path(results_directory, "fit_diagnostics.csv")
)

# Step 3: run independent checks on the data, timing, likelihood, and scores.
validation <- tryCatch(
  validate_analysis(
    merged_data = merged_data,
    predictions = predictions,
    diagnostics = fit_diagnostics,
    metrics = metrics,
    input_path = input_file,
    config = config
  ),
  error = function(error) error
)

if (inherits(validation, "error")) {
  writeLines(
    conditionMessage(validation),
    file.path(results_directory, "validation_failure.txt")
  )
  stop(validation)
}

validation_file <- file.path(results_directory, "validation_failure.txt")
if (file.exists(validation_file)) unlink(validation_file)
write_project_csv(
  validation,
  file.path(results_directory, "validation_checks.csv")
)

# Step 4: compare all prespecified reduced predictor sets on the same forecast
# dates and rolling training windows as the full five-predictor model.
reduced_comparison <- run_reduced_model_comparison(
  merged_data = merged_data,
  config = config,
  baseline_predictions = predictions,
  baseline_diagnostics = fit_diagnostics,
  output_dir = results_directory
)

# Save session details and all main objects for reproducibility.
session_information <- utils::sessionInfo()
writeLines(
  capture.output(session_information),
  file.path(results_directory, "session_info.txt")
)
saveRDS(
  list(
    config = config,
    predictions = predictions,
    metrics = metrics,
    diagnostics = fit_diagnostics,
    validation = validation,
    reduced_model_comparison = reduced_comparison,
    session_info = session_information
  ),
  file.path(results_directory, "analysis_objects.rds")
)

# Step 5: run collinearity diagnostics, presentation charts/tables, and the
# additional figures and tables used in the written report.
run_collinearity_diagnostics(
  merged_data = merged_data,
  predictions = predictions,
  config = config,
  output_dir = results_directory
)

create_presentation_materials(output_dir = output_directory)
create_report_materials(output_dir = output_directory)

message(" ")
message("Analysis completed.")
message("Core results:          ", results_directory)
message("Presentation materials: ", file.path(output_directory, "presentation_materials"))
message("Report figures:         ", file.path(output_directory, "report_figures"))
message(sprintf(
  "Exceedances: %d | Rolling prediction rows: %d | Validation checks passed: %d/%d",
  sum(merged_data$exceedance == 1L, na.rm = TRUE),
  nrow(predictions),
  sum(validation$passed),
  nrow(validation)
))
message(sprintf(
  "Reduced-model forecasts: %d | Comparison checks passed: %d/%d",
  nrow(reduced_comparison$predictions),
  sum(reduced_comparison$validation$passed),
  nrow(reduced_comparison$validation)
))
