# Create presentation-ready tables and charts from the saved model outputs.

create_presentation_materials <- function(output_dir = ".") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install ggplot2 before creating presentation materials.")
  }

  results_dir <- file.path(output_dir, "results")
  presentation_dir <- file.path(output_dir, "presentation_materials")
  dir.create(presentation_dir, showWarnings = FALSE, recursive = TRUE)

  read_result <- function(filename) {
    path <- file.path(results_dir, filename)
    if (!file.exists(path)) stop("Required result file not found: ", path)
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }

  nice_names <- c(
    constant_gpd = "Constant GPD baseline (no covariates)",
    full_5 = paste(
      "Nifty volatility + S&P 500 volatility + Gold volatility +",
      "India VIX + USD/INR volatility"
    ),
    A = "Nifty volatility + S&P 500 volatility + Gold volatility + USD/INR volatility",
    B = "S&P 500 volatility + Gold volatility + India VIX + USD/INR volatility",
    C = "Nifty 50 volatility + USD/INR volatility",
    D = "India VIX + USD/INR volatility"
  )

  predictor_names <- c(
    "nifty_vol_lag1", "sp500_vol_lag1", "gold_vol_lag1",
    "india_vix_lag1", "usdinr_vol_lag1"
  )
  predictor_labels <- c(
    nifty_vol_lag1 = "Nifty volatility",
    sp500_vol_lag1 = "S&P 500 volatility",
    gold_vol_lag1 = "Gold volatility",
    india_vix_lag1 = "India VIX",
    usdinr_vol_lag1 = "USD/INR volatility"
  )

  presentation_theme <- function(base_size = 12) {
    ggplot2::theme_minimal(base_size = base_size) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.2)),
        plot.subtitle = ggplot2::element_text(colour = "#4B5563"),
        axis.title = ggplot2::element_text(face = "bold"),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "bottom",
        legend.title = ggplot2::element_blank()
      )
  }

  save_plot <- function(plot, filename, width = 11, height = 6) {
    ggplot2::ggsave(
      filename = file.path(presentation_dir, filename),
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = 300,
      bg = "white"
    )
  }

  merged <- read_result("merged_data.csv")
  merged$Date <- as.Date(merged$Date)
  analysis_objects <- readRDS(file.path(results_dir, "analysis_objects.rds"))
  rolling_window_size <- analysis_objects$config$rolling_observations
  event_count <- sum(as.numeric(merged$exceedance) == 1, na.rm = TRUE)
  threshold_values <- unique(stats::na.omit(as.numeric(merged$threshold)))
  threshold_value <- if (length(threshold_values) == 1L) threshold_values else 2

  comparison_predictions <- read_result("reduced_model_comparison_predictions.csv")
  comparison_predictions$Date <- as.Date(comparison_predictions$Date)
  forecast_dates <- sort(unique(comparison_predictions$Date))

  # Table 1: data coverage and rolling evaluation design.
  data_summary <- data.frame(
    Item = c(
      "Nifty 50 data start",
      "Nifty 50 data end",
      "Extreme-loss threshold (percent)",
      "Extreme-loss exceedances",
      "Rolling forecast dates",
      "Prior Nifty trading dates per rolling window",
      "Forecast target",
      "Forecast scope"
    ),
    Value = c(
      format(min(merged$Date, na.rm = TRUE)),
      format(max(merged$Date, na.rm = TRUE)),
      paste(threshold_value, collapse = ", "),
      event_count,
      length(forecast_dates),
      rolling_window_size,
      "Excess loss above the threshold",
      "Severity conditional on threshold exceedance"
    ),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    data_summary,
    file.path(presentation_dir, "table_data_and_evaluation_summary.csv"),
    row.names = FALSE
  )

  # Calculate forecast R-squared directly for each model-set/prior group.
  grouped_predictions <- split(
    comparison_predictions,
    interaction(
      comparison_predictions$model_set,
      comparison_predictions$prior,
      drop = TRUE
    )
  )
  forecast_r2 <- do.call(rbind, lapply(grouped_predictions, function(data) {
    observed <- as.numeric(data$observed_excess)
    predicted <- as.numeric(data$predicted_mean_excess)
    denominator <- sum((observed - mean(observed))^2)

    data.frame(
      model_set = data$model_set[1],
      prior = data$prior[1],
      forecast_R_squared = if (denominator > 0) {
        1 - sum((observed - predicted)^2) / denominator
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }))

  # Table 2: performance measures and rolling VIF summaries.
  metrics <- read_result("reduced_model_comparison_metrics.csv")
  fit_diagnostics <- read_result("reduced_model_comparison_fit_diagnostics.csv")

  vif_rows <- unique(
    fit_diagnostics[c("forecast_date", "model_set", "max_vif")]
  )
  vif_rows$max_vif <- as.numeric(vif_rows$max_vif)
  vif_rows <- vif_rows[
    vif_rows$model_set != "constant_gpd" & is.finite(vif_rows$max_vif),
    , drop = FALSE
  ]
  vif_summary <- if (nrow(vif_rows) > 0L) {
    do.call(rbind, lapply(split(vif_rows, vif_rows$model_set), function(data) {
      data.frame(
        model_set = data$model_set[1],
        median_rolling_max_VIF = stats::median(data$max_vif),
        maximum_rolling_VIF = max(data$max_vif),
        stringsAsFactors = FALSE
      )
    }))
  } else {
    data.frame(
      model_set = character(),
      median_rolling_max_VIF = numeric(),
      maximum_rolling_VIF = numeric()
    )
  }

  performance <- merge(
    metrics, forecast_r2,
    by = c("model_set", "prior"), all.x = TRUE, sort = FALSE
  )
  performance <- merge(
    performance, vif_summary,
    by = "model_set", all.x = TRUE, sort = FALSE
  )
  performance$Model <- unname(nice_names[performance$model_set])
  performance$Prior <- performance$prior
  performance$Forecasts <- performance$forecasts

  presentation_metrics <- performance[c(
    "Model", "Prior", "Forecasts", "RMSE", "MAE",
    "forecast_R_squared", "mean_conditional_log_score",
    "median_rolling_max_VIF", "maximum_rolling_VIF"
  )]
  names(presentation_metrics) <- c(
    "Model", "Prior", "Forecast dates", "RMSE", "MAE",
    "Out-of-sample R-squared", "Mean conditional log score",
    "Median rolling max VIF", "Maximum rolling VIF"
  )
  utils::write.csv(
    presentation_metrics,
    file.path(presentation_dir, "table_all_model_results.csv"),
    row.names = FALSE
  )

  # Table 3: predictor definitions.
  model_definitions <- data.frame(
    Model = unname(nice_names[c("constant_gpd", "full_5", "A", "B", "C", "D")]),
    Covariates = c(
      "None; intercept only",
      "Nifty volatility; S&P 500 volatility; Gold volatility; India VIX; USD/INR volatility",
      "Nifty volatility; S&P 500 volatility; Gold volatility; USD/INR volatility",
      "S&P 500 volatility; Gold volatility; India VIX; USD/INR volatility",
      "Nifty volatility; USD/INR volatility",
      "India VIX; USD/INR volatility"
    ),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    model_definitions,
    file.path(presentation_dir, "table_model_definitions.csv"),
    row.names = FALSE
  )

  # Table 4: maximum VIF by model set.
  vif_summary$Model <- unname(nice_names[vif_summary$model_set])
  vif_table <- vif_summary[c(
    "Model", "median_rolling_max_VIF", "maximum_rolling_VIF"
  )]
  names(vif_table) <- c(
    "Model", "Median rolling max VIF", "Maximum rolling VIF"
  )
  utils::write.csv(
    vif_table,
    file.path(presentation_dir, "table_rolling_VIF_summary.csv"),
    row.names = FALSE
  )

  # Figure 1: Nifty daily returns and the extreme-loss threshold.
  return_data <- merged[
    !is.na(merged$nifty_return),
    c("Date", "nifty_return", "exceedance")
  ]
  return_data$exceedance <- as.numeric(return_data$exceedance)

  plot_returns <- ggplot2::ggplot(
    return_data,
    ggplot2::aes(x = Date, y = nifty_return)
  ) +
    ggplot2::geom_hline(
      yintercept = -threshold_value,
      colour = "#C0392B", linetype = "dashed", linewidth = 0.8
    ) +
    ggplot2::geom_line(
      colour = "#526D82", linewidth = 0.25, alpha = 0.8
    ) +
    ggplot2::geom_point(
      data = return_data[return_data$exceedance == 1, ],
      colour = "#C0392B", size = 1.1, alpha = 0.75
    ) +
    ggplot2::labs(
      title = "Nifty 50 daily returns and extreme-loss threshold",
      subtitle = paste0(
        event_count, " threshold exceedances; dashed line marks a loss of ",
        threshold_value, "%"
      ),
      x = NULL,
      y = "Daily return (%)"
    ) +
    presentation_theme()
  save_plot(plot_returns, "graph_1_nifty_returns_threshold.png", 12, 5.5)

  # Figure 2: observed severity and selected model forecasts on common dates.
  chosen <- comparison_predictions[
    (comparison_predictions$model_set == "constant_gpd" &
       comparison_predictions$prior == "baseline") |
      (comparison_predictions$model_set == "C" &
         comparison_predictions$prior == "g-prior") |
      (comparison_predictions$model_set == "D" &
         comparison_predictions$prior == "cauchy"),
    , drop = FALSE
  ]
  chosen <- chosen[order(chosen$Date, chosen$model_set), ]

  observed_line <- unique(chosen[c("Date", "observed_excess")])
  observed_line$Series <- "Observed excess loss"
  observed_line$Value <- observed_line$observed_excess

  predicted_lines <- chosen[c(
    "Date", "model_set", "prior", "predicted_mean_excess"
  )]
  predicted_lines$Series <- ifelse(
    predicted_lines$model_set == "constant_gpd",
    "Constant GPD baseline",
    ifelse(
      predicted_lines$model_set == "C",
      "Nifty volatility + USD/INR volatility (g-prior)",
      "India VIX + USD/INR volatility (Cauchy)"
    )
  )
  predicted_lines$Value <- predicted_lines$predicted_mean_excess

  plot_lines <- rbind(
    observed_line[c("Date", "Series", "Value")],
    predicted_lines[c("Date", "Series", "Value")]
  )
  plot_forecasts <- ggplot2::ggplot(
    plot_lines,
    ggplot2::aes(x = Date, y = Value, colour = Series, group = Series)
  ) +
    ggplot2::geom_line(linewidth = 0.75, alpha = 0.9) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_colour_manual(values = c(
      "Observed excess loss" = "#111827",
      "Constant GPD baseline" = "#9CA3AF",
      "Nifty volatility + USD/INR volatility (g-prior)" = "#2563EB",
      "India VIX + USD/INR volatility (Cauchy)" = "#D97706"
    )) +
    ggplot2::labs(
      title = "Observed and predicted extreme-loss severity",
      subtitle = paste(
        "Conditional mean excess forecasts on the same",
        length(forecast_dates), "rolling forecast dates"
      ),
      x = NULL,
      y = "Excess loss above threshold (percentage points)"
    ) +
    presentation_theme(11) +
    ggplot2::theme(legend.text = ggplot2::element_text(size = 9))
  save_plot(plot_forecasts, "graph_2_actual_vs_predicted_severity.png", 12, 6)

  # Figure 3: distribution of fold-level maximum VIF values.
  vif_rows$Model <- unname(nice_names[vif_rows$model_set])
  model_order <- vif_table$Model[order(vif_table[["Maximum rolling max VIF"]])]
  vif_rows$Model <- factor(vif_rows$Model, levels = model_order)

  plot_vif <- ggplot2::ggplot(vif_rows, ggplot2::aes(x = max_vif, y = Model)) +
    ggplot2::geom_vline(
      xintercept = 5, linetype = "dashed", colour = "#D97706", linewidth = 0.7
    ) +
    ggplot2::geom_vline(
      xintercept = 10, linetype = "dotted", colour = "#C0392B", linewidth = 0.7
    ) +
    ggplot2::geom_boxplot(
      fill = "#8FB9A8", colour = "#355C50",
      outlier.colour = "#C0392B", outlier.size = 1.1
    ) +
    ggplot2::labs(
      title = "Rolling multicollinearity by predictor set",
      subtitle = "Each box summarises fold-specific maximum VIF; 5 and 10 are common heuristics",
      x = "Maximum VIF within training fold",
      y = NULL
    ) +
    presentation_theme(10)
  save_plot(plot_vif, "graph_3_rolling_VIF_by_model.png", 12, 6.5)

  # Figure 4: predictor correlation heatmap on complete exceedance rows.
  complete_x <- merged[
    as.numeric(merged$exceedance) == 1 &
      is.finite(merged$excess) &
      merged$excess > 0 &
      stats::complete.cases(merged[, predictor_names]),
    predictor_names,
    drop = FALSE
  ]
  correlation <- stats::cor(complete_x)
  utils::write.csv(
    correlation,
    file.path(presentation_dir, "table_predictor_correlation_matrix.csv"),
    row.names = TRUE
  )

  correlation_long <- as.data.frame(
    as.table(correlation), stringsAsFactors = FALSE
  )
  names(correlation_long) <- c("Predictor_1", "Predictor_2", "Correlation")
  correlation_long$Predictor_1 <- unname(
    predictor_labels[as.character(correlation_long$Predictor_1)]
  )
  correlation_long$Predictor_2 <- unname(
    predictor_labels[as.character(correlation_long$Predictor_2)]
  )

  plot_correlation <- ggplot2::ggplot(
    correlation_long,
    ggplot2::aes(x = Predictor_1, y = Predictor_2, fill = Correlation)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", Correlation)),
      size = 3.5, colour = "#111827"
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-1, 1)
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "Correlation among lagged predictors",
      subtitle = paste(
        "Pearson correlations on", nrow(complete_x),
        "exceedances complete for all five covariates"
      ),
      x = NULL, y = NULL, fill = "Correlation"
    ) +
    presentation_theme(10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  save_plot(plot_correlation, "graph_4_predictor_correlation_heatmap.png", 9, 7)

  # Figure 5: error comparison for selected specifications.
  shortlist <- performance[
    (performance$model_set == "constant_gpd" & performance$prior == "baseline") |
      (performance$model_set == "full_5" & performance$prior == "g-prior") |
      (performance$model_set == "C" & performance$prior == "g-prior") |
      (performance$model_set == "D" & performance$prior == "cauchy"),
    , drop = FALSE
  ]
  shortlist$Specification <- ifelse(
    shortlist$model_set == "constant_gpd",
    "Constant GPD baseline",
    ifelse(
      shortlist$model_set == "full_5",
      "Full five-covariate model (g-prior)",
      ifelse(
        shortlist$model_set == "C",
        "Nifty volatility + USD/INR volatility (g-prior)",
        "India VIX + USD/INR volatility (Cauchy)"
      )
    )
  )

  error_data <- rbind(
    data.frame(
      Specification = shortlist$Specification,
      Metric = "RMSE",
      Value = as.numeric(shortlist$RMSE)
    ),
    data.frame(
      Specification = shortlist$Specification,
      Metric = "MAE",
      Value = as.numeric(shortlist$MAE)
    )
  )
  error_data$Specification <- factor(
    error_data$Specification,
    levels = rev(shortlist$Specification)
  )

  plot_errors <- ggplot2::ggplot(
    error_data,
    ggplot2::aes(x = Value, y = Specification, fill = Specification)
  ) +
    ggplot2::geom_col(show.legend = FALSE, width = 0.68) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.3f", Value)),
      hjust = -0.15, size = 3.4
    ) +
    ggplot2::facet_wrap(~Metric, scales = "free_x") +
    ggplot2::scale_fill_manual(values = c(
      "Constant GPD baseline" = "#9CA3AF",
      "Full five-covariate model (g-prior)" = "#566573",
      "Nifty volatility + USD/INR volatility (g-prior)" = "#2563EB",
      "India VIX + USD/INR volatility (Cauchy)" = "#D97706"
    )) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.14))
    ) +
    ggplot2::labs(
      title = "Rolling forecast error for selected specifications",
      subtitle = paste(
        "All values use the same", length(forecast_dates),
        "forecast dates; lower values indicate smaller errors"
      ),
      x = "Forecast error",
      y = NULL
    ) +
    presentation_theme(10) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
  save_plot(plot_errors, "graph_5_selected_model_errors.png", 13, 6.5)

  notes <- c(
    "# Presentation materials",
    "",
    "Generated from the project's saved rolling forecasts and diagnostics by `R/presentation_materials.R`.",
    "",
    "## Files",
    "- `graph_1_nifty_returns_threshold.png`: Nifty 50 daily returns, threshold, and exceedance dates.",
    "- `graph_2_actual_vs_predicted_severity.png`: observed excess severity and selected conditional-mean forecasts.",
    "- `graph_3_rolling_VIF_by_model.png`: distribution of fold-specific maximum VIF by predictor set.",
    "- `graph_4_predictor_correlation_heatmap.png`: Pearson correlation heatmap for the five lagged predictors.",
    "- `graph_5_selected_model_errors.png`: RMSE and MAE comparison for selected model specifications.",
    "- CSV tables: data summary, model results, model definitions, rolling VIF summary, and correlation matrix.",
    "",
    "## Metric definitions and cautions",
    "- RMSE and MAE compare predicted conditional mean excess loss with observed excess loss on held-out exceedance dates.",
    "- Out-of-sample R-squared is 1 - SSE/SST, with SST based on the mean observed excess across the forecast dates. It is not adjusted R-squared and can be negative.",
    "- A larger (less negative) conditional log score is better. It assesses the fitted excess-loss distribution conditional on an exceedance.",
    "- VIF diagnoses predictor collinearity. It does not, by itself, establish predictive overfitting.",
    "- The model predicts excess-loss severity conditional on a threshold exceedance, not the probability that an exceedance occurs.",
    "- The same forecast dates are used for model-set comparisons. Treat a selected winner as exploratory unless a separate untouched evaluation period is retained."
  )
  writeLines(notes, file.path(presentation_dir, "README.md"))

  message("Presentation materials saved to: ", normalizePath(presentation_dir))
  invisible(list(
    data_summary = data_summary,
    metrics = presentation_metrics,
    model_definitions = model_definitions,
    vif_summary = vif_table
  ))
}
