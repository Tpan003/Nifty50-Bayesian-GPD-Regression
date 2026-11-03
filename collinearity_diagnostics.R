# Full-sample and rolling-window multicollinearity diagnostics for the five
# lagged predictors used in the GPD severity models.

run_collinearity_diagnostics <- function(merged_data, predictions, config,
                                         output_dir = "results") {
  predictors <- c(
    "nifty_vol_lag1", "sp500_vol_lag1", "gold_vol_lag1",
    "india_vix_lag1", "usdinr_vol_lag1"
  )
  required <- c("Date", "exceedance", "excess", predictors)
  if (!all(required %in% names(merged_data))) {
    stop("Merged data is missing columns required for collinearity diagnostics.")
  }
  if (nrow(predictions) == 0L) {
    stop("No rolling forecast dates are available for collinearity diagnostics.")
  }

  merged_data$Date <- as.Date(merged_data$Date)
  predictions$Date <- as.Date(predictions$Date)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  profile_design <- function(data) {
    x <- as.data.frame(lapply(data[, predictors, drop = FALSE], as.numeric))
    standard_deviations <- vapply(x, stats::sd, numeric(1))

    if (nrow(x) < 2L || any(!is.finite(standard_deviations)) ||
        any(standard_deviations <= 0)) {
      stop("A collinearity design has too few rows or a zero-variance predictor.")
    }

    correlation <- stats::cor(x, method = "pearson")
    upper <- which(upper.tri(correlation), arr.ind = TRUE)
    pair_values <- correlation[upper]
    strongest_pair <- which.max(abs(pair_values))

    vif <- vapply(predictors, function(response) {
      other_predictors <- setdiff(predictors, response)
      fit <- stats::lm(
        stats::reformulate(other_predictors, response = response),
        data = x
      )
      r_squared <- summary(fit)$r.squared
      if (!is.finite(r_squared) || r_squared >= 1) {
        Inf
      } else {
        1 / (1 - r_squared)
      }
    }, numeric(1))
    names(vif) <- paste0("vif_", predictors)

    eigenvalues <- eigen(
      correlation, symmetric = TRUE, only.values = TRUE
    )$values
    smallest_eigenvalue <- min(eigenvalues)
    condition_index <- if (smallest_eigenvalue <= 0) {
      Inf
    } else {
      sqrt(max(eigenvalues) / smallest_eigenvalue)
    }

    list(
      summary = data.frame(
        training_exceedances = nrow(x),
        design_rank = qr(cbind(Intercept = 1, as.matrix(x)))$rank,
        max_abs_pairwise_r = abs(pair_values[strongest_pair]),
        max_correlated_pair = paste(
          rownames(correlation)[upper[strongest_pair, 1]],
          colnames(correlation)[upper[strongest_pair, 2]],
          sep = " -- "
        ),
        max_vif = max(vif),
        max_vif_predictor = sub(
          "^vif_", "", names(vif)[which.max(vif)]
        ),
        condition_index = condition_index,
        stringsAsFactors = FALSE
      ),
      correlation = correlation,
      vif = vif
    )
  }

  # Descriptive full-sample diagnostic on complete exceedance rows.
  full_events <- merged_data[
    merged_data$exceedance == 1L &
      is.finite(merged_data$excess) &
      merged_data$excess > 0 &
      complete.cases(merged_data[, predictors, drop = FALSE]),
    , drop = FALSE
  ]
  full_profile <- profile_design(full_events)

  utils::write.csv(
    full_profile$correlation,
    file.path(output_dir, "collinearity_full_sample_pearson.csv"),
    row.names = TRUE
  )
  utils::write.csv(
    stats::cor(full_events[, predictors, drop = FALSE], method = "spearman"),
    file.path(output_dir, "collinearity_full_sample_spearman.csv"),
    row.names = TRUE
  )

  full_vif <- data.frame(
    predictor = predictors,
    VIF = as.numeric(full_profile$vif),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    full_vif,
    file.path(output_dir, "collinearity_full_sample_vif.csv"),
    row.names = FALSE
  )

  # Recreate the exact preceding-date training window for every scored date.
  forecast_dates <- sort(unique(predictions$Date))
  rolling_rows <- vector("list", length(forecast_dates))
  pair_rows <- vector("list", length(forecast_dates))

  for (i in seq_along(forecast_dates)) {
    forecast_date <- forecast_dates[i]
    previous_dates <- sort(unique(
      merged_data$Date[merged_data$Date < forecast_date]
    ))
    window_dates <- tail(previous_dates, config$rolling_observations)

    if (length(window_dates) != config$rolling_observations ||
        max(window_dates) >= forecast_date) {
      stop("A rolling collinearity window does not match the model window.")
    }

    window_data <- merged_data[
      merged_data$Date %in% window_dates, , drop = FALSE
    ]
    training_events <- window_data[
      window_data$exceedance == 1L &
        is.finite(window_data$excess) &
        window_data$excess > 0 &
        complete.cases(window_data[, predictors, drop = FALSE]),
      , drop = FALSE
    ]

    if (nrow(training_events) < config$minimum_training_exceedances) {
      stop("A rolling collinearity window has too few training exceedances.")
    }

    profile <- profile_design(training_events)
    rolling_rows[[i]] <- cbind(
      data.frame(
        forecast_date = forecast_date,
        window_first_date = min(window_dates),
        window_last_date = max(window_dates),
        stringsAsFactors = FALSE
      ),
      profile$summary,
      as.data.frame(as.list(profile$vif), check.names = FALSE)
    )

    upper <- which(upper.tri(profile$correlation), arr.ind = TRUE)
    pair_rows[[i]] <- data.frame(
      forecast_date = forecast_date,
      predictor_1 = rownames(profile$correlation)[upper[, 1]],
      predictor_2 = colnames(profile$correlation)[upper[, 2]],
      pearson_r = profile$correlation[upper],
      stringsAsFactors = FALSE
    )
  }

  rolling <- do.call(rbind, rolling_rows)
  pairs <- do.call(rbind, pair_rows)

  if (any(rolling$design_rank != length(predictors) + 1L) ||
      any(rolling$training_exceedances < config$minimum_training_exceedances) ||
      any(rolling$window_last_date >= rolling$forecast_date)) {
    stop("A rolling collinearity check failed rank, sample-size, or date-order validation.")
  }

  utils::write.csv(
    rolling,
    file.path(output_dir, "collinearity_rolling_windows.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    pairs,
    file.path(output_dir, "collinearity_rolling_pairwise.csv"),
    row.names = FALSE
  )

  pair_summary <- do.call(rbind, lapply(
    split(pairs, interaction(pairs$predictor_1, pairs$predictor_2, drop = TRUE)),
    function(data) {
      first <- data$predictor_1[1]
      second <- data$predictor_2[1]
      data.frame(
        predictor_1 = first,
        predictor_2 = second,
        full_sample_pearson_r = full_profile$correlation[first, second],
        rolling_median_r = stats::median(data$pearson_r),
        rolling_min_r = min(data$pearson_r),
        rolling_max_r = max(data$pearson_r),
        rolling_max_abs_r = max(abs(data$pearson_r)),
        stringsAsFactors = FALSE
      )
    }
  ))
  pair_summary <- pair_summary[order(-pair_summary$rolling_max_abs_r), ]
  utils::write.csv(
    pair_summary,
    file.path(output_dir, "collinearity_pair_summary.csv"),
    row.names = FALSE
  )

  summary_row <- data.frame(
    full_sample_complete_exceedances = nrow(full_events),
    rolling_forecast_dates = nrow(rolling),
    min_rolling_training_exceedances = min(rolling$training_exceedances),
    max_rolling_training_exceedances = max(rolling$training_exceedances),
    median_rolling_max_abs_r = stats::median(rolling$max_abs_pairwise_r),
    maximum_rolling_max_abs_r = max(rolling$max_abs_pairwise_r),
    windows_max_abs_r_ge_0_8 = sum(rolling$max_abs_pairwise_r >= 0.8),
    median_rolling_max_vif = stats::median(rolling$max_vif),
    maximum_rolling_vif = max(rolling$max_vif),
    windows_max_vif_ge_5 = sum(rolling$max_vif >= 5),
    windows_max_vif_ge_10 = sum(rolling$max_vif >= 10),
    median_rolling_condition_index = stats::median(rolling$condition_index),
    maximum_rolling_condition_index = max(rolling$condition_index),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    summary_row,
    file.path(output_dir, "collinearity_summary.csv"),
    row.names = FALSE
  )

  rolling_vif <- do.call(rbind, lapply(predictors, function(variable) {
    values <- rolling[[paste0("vif_", variable)]]
    data.frame(
      predictor = variable,
      median_vif = stats::median(values),
      maximum_vif = max(values),
      windows_vif_ge_5 = sum(values >= 5),
      windows_vif_ge_10 = sum(values >= 10),
      stringsAsFactors = FALSE
    )
  }))

  full_upper <- which(upper.tri(full_profile$correlation), arr.ind = TRUE)
  full_pairs <- data.frame(
    predictor_1 = rownames(full_profile$correlation)[full_upper[, 1]],
    predictor_2 = colnames(full_profile$correlation)[full_upper[, 2]],
    pearson_r = full_profile$correlation[full_upper],
    stringsAsFactors = FALSE
  )
  full_pairs <- full_pairs[order(-abs(full_pairs$pearson_r)), ]
  top_full_pairs <- utils::head(full_pairs, 5)
  top_rolling_pairs <- utils::head(pair_summary, 5)

  vif_table_lines <- c(
    "| Predictor | Full-sample VIF | Median rolling VIF | Maximum rolling VIF | Folds with VIF >= 5 | Folds with VIF >= 10 |",
    "|:--|--:|--:|--:|--:|--:|",
    vapply(seq_len(nrow(rolling_vif)), function(i) {
      variable <- rolling_vif$predictor[i]
      full_value <- full_vif$VIF[match(variable, full_vif$predictor)]
      sprintf(
        "| %s | %.2f | %.2f | %.2f | %d/%d | %d/%d |",
        variable, full_value, rolling_vif$median_vif[i],
        rolling_vif$maximum_vif[i], rolling_vif$windows_vif_ge_5[i],
        nrow(rolling), rolling_vif$windows_vif_ge_10[i], nrow(rolling)
      )
    }, character(1))
  )

  full_pair_lines <- c(
    "| Predictor pair | Pearson r |",
    "|:--|--:|",
    vapply(seq_len(nrow(top_full_pairs)), function(i) {
      sprintf(
        "| %s / %s | %.3f |",
        top_full_pairs$predictor_1[i], top_full_pairs$predictor_2[i],
        top_full_pairs$pearson_r[i]
      )
    }, character(1))
  )

  rolling_pair_lines <- c(
    "| Predictor pair | Median rolling r | Maximum rolling absolute r |",
    "|:--|--:|--:|",
    vapply(seq_len(nrow(top_rolling_pairs)), function(i) {
      sprintf(
        "| %s / %s | %.3f | %.3f |",
        top_rolling_pairs$predictor_1[i], top_rolling_pairs$predictor_2[i],
        top_rolling_pairs$rolling_median_r[i],
        top_rolling_pairs$rolling_max_abs_r[i]
      )
    }, character(1))
  )

  # Summarise model scores for context. These are conditional severity scores,
  # not scores for the probability of an exceedance occurring.
  model_metrics_path <- file.path(output_dir, "rolling_metrics.csv")
  score_note <- "The rolling metrics file was not found, so predictive scores were not added to this report."
  if (file.exists(model_metrics_path)) {
    model_metrics <- utils::read.csv(model_metrics_path)
    regression_metrics <- model_metrics[model_metrics$model != "constant_gpd", ]
    baseline_metrics <- model_metrics[model_metrics$model == "constant_gpd", ]
    if (nrow(regression_metrics) > 0L && nrow(baseline_metrics) == 1L) {
      score_note <- sprintf(
        paste0(
          "Across %d forecast dates, the covariate models have RMSE values from %.3f to %.3f, ",
          "compared with %.3f for the constant GPD baseline. Their mean conditional log scores ",
          "range from %.3f to %.3f, compared with %.3f for the baseline. Higher log scores are better."
        ),
        nrow(rolling), min(regression_metrics$RMSE, na.rm = TRUE),
        max(regression_metrics$RMSE, na.rm = TRUE), baseline_metrics$RMSE,
        min(regression_metrics$mean_conditional_log_score, na.rm = TRUE),
        max(regression_metrics$mean_conditional_log_score, na.rm = TRUE),
        baseline_metrics$mean_conditional_log_score
      )
    }
  }

  report_lines <- c(
    "# Predictor multicollinearity diagnostics",
    "",
    "This report checks the five lagged predictors on complete Nifty 50 exceedance observations. The rolling analysis reconstructs the preceding-date training window for each held-out exceedance date.",
    "",
    sprintf(
      "The full-sample check uses %d complete exceedances. The rolling audit covers %d forecast dates, with %d to %d training exceedances per fold and %d prior Nifty trading dates in each window.",
      nrow(full_events), nrow(rolling), min(rolling$training_exceedances),
      max(rolling$training_exceedances), config$rolling_observations
    ),
    "",
    "## Main findings",
    "",
    sprintf(
      "- The strongest full-sample Pearson association is **%s / %s (r = %.3f)**.",
      top_full_pairs$predictor_1[1], top_full_pairs$predictor_2[1],
      top_full_pairs$pearson_r[1]
    ),
    sprintf(
      "- The median of the fold-level maximum absolute correlations is **%.3f**; the largest is **%.3f**. %d of %d folds have a maximum absolute correlation of at least 0.8.",
      summary_row$median_rolling_max_abs_r,
      summary_row$maximum_rolling_max_abs_r,
      summary_row$windows_max_abs_r_ge_0_8, nrow(rolling)
    ),
    sprintf(
      "- The median fold maximum VIF is **%.2f**, and the maximum is **%.2f**. %d folds have a maximum VIF of at least 5; %d folds have a maximum VIF of at least 10.",
      summary_row$median_rolling_max_vif,
      summary_row$maximum_rolling_max_vif,
      summary_row$windows_max_vif_ge_5,
      summary_row$windows_max_vif_ge_10
    ),
    sprintf(
      "- Every rolling design has rank %d (an intercept and five predictors). The maximum condition index is **%.2f**.",
      length(predictors) + 1L, summary_row$maximum_rolling_condition_index
    ),
    "",
    "VIF values around 5 and 10 are commonly used as warning thresholds, not formal pass/fail tests. Correlation and VIF are unchanged by linear rescaling, so the raw lagged predictors give the same collinearity diagnostics as standardised predictors.",
    "",
    "## VIF by predictor",
    "",
    vif_table_lines,
    "",
    "## Strongest full-sample correlations",
    "",
    full_pair_lines,
    "",
    "## Strongest relationships across rolling folds",
    "",
    rolling_pair_lines,
    "",
    "## Interpretation",
    "",
    "High correlation and VIF indicate overlapping information and a risk of unstable regression coefficients. They do not, by themselves, prove that forecasts are overfit. The rolling training samples contain a limited number of exceedances relative to the number of parameters, so regularisation and reduced-model comparisons are useful sensitivity checks.",
    "",
    score_note,
    "",
    "The log scores are conditional on a threshold exceedance. This model evaluates excess-loss severity after an exceedance; it does not model the probability that an exceedance occurs.",
    "",
    "## Recommended follow-up",
    "",
    "1. Compare the prespecified reduced predictor sets on the same rolling forecast dates. Do not remove a variable based only on its correlation with another variable.",
    "2. Use RMSE and MAE for conditional-mean accuracy and the conditional log score for the fitted tail distribution. Consider PIT diagnostics as a separate calibration check.",
    "3. If predictor selection or prior tuning is based on rolling results, use nested rolling validation and retain a final untouched evaluation period.",
    "",
    "## Reproduce",
    "",
    "Run `Rscript run_analysis.R <workbook.xlsx> <output-folder>` from the repository root. The diagnostics are generated automatically and saved as `results/collinearity_*.csv` and `results/collinearity_report.md`."
  )

  writeLines(
    report_lines,
    file.path(output_dir, "collinearity_report.md"),
    useBytes = TRUE
  )

  invisible(list(
    full_sample = full_profile,
    rolling = rolling,
    pairwise = pairs,
    pair_summary = pair_summary,
    summary = summary_row,
    vif = full_vif
  ))
}
