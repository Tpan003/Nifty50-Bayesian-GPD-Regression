# Build the additional report tables and figures used in the written report.
# These calculations use base R graphics so that no plotting package is needed.

create_report_materials <- function(output_dir = ".") {
  results_dir <- file.path(output_dir, "results")
  figure_dir <- file.path(output_dir, "report_figures")
  dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

  read_result <- function(filename) {
    path <- file.path(results_dir, filename)
    if (!file.exists(path)) stop("Required result file not found: ", path)
    utils::read.csv(path, stringsAsFactors = FALSE)
  }

  write_result <- function(data, filename) {
    utils::write.csv(
      data,
      file.path(results_dir, filename),
      row.names = FALSE,
      na = ""
    )
  }

  open_figure <- function(filename, width = 10, height = 5.5) {
    grDevices::png(
      file.path(figure_dir, filename),
      width = width,
      height = height,
      units = "in",
      res = 300,
      bg = "white"
    )
  }

  config <- project_config()
  prior_spec <- config$prior_specification
  threshold <- config$threshold
  model_sets <- reduced_model_sets()
  full_predictors <- model_sets$full_5

  merged <- read_result("merged_data.csv")
  merged$Date <- as.Date(merged$Date)
  exceedances <- merged[
    merged$exceedance == 1 &
      is.finite(merged$excess) &
      merged$excess > 0,
    , drop = FALSE
  ]
  complete_exceedances <- exceedances[
    stats::complete.cases(exceedances[, full_predictors, drop = FALSE]),
    , drop = FALSE
  ]

  comparison_predictions <- read_result("reduced_model_comparison_predictions.csv")
  comparison_predictions$Date <- as.Date(comparison_predictions$Date)
  comparison_metrics <- read_result("reduced_model_comparison_metrics.csv")
  comparison_windows <- read_result("reduced_model_comparison_windows.csv")
  comparison_windows$forecast_date <- as.Date(comparison_windows$forecast_date)

  excess <- exceedances$excess
  if (length(excess) < 2L) {
    stop("There are too few exceedances to create the report tables and figures.")
  }

  # GPD utilities used for distributional summaries and diagnostics.
  fit_pooled_gpd_mle <- function(z) {
    negative_log_likelihood <- function(parameters) {
      value <- -sum(gpd_log_density(
        z, rep(parameters[1], length(z)), parameters[2]
      ))
      if (is.finite(value)) value else 1e100
    }

    fit <- stats::nlminb(
      start = c(log(mean(z)), 0.1),
      objective = negative_log_likelihood,
      lower = c(-10, -5),
      upper = c(10, 1 - 1e-6)
    )
    if (fit$convergence != 0L || !is.finite(fit$objective) || fit$objective >= 1e99) {
      stop("The pooled GPD maximum-likelihood fit did not converge.")
    }

    list(
      sigma = exp(fit$par[1]),
      xi = fit$par[2],
      loglik = -fit$objective,
      convergence = fit$convergence
    )
  }

  gpd_quantile <- function(probability, scale, shape) {
    if (abs(shape) < 1e-7) {
      -scale * log(1 - probability)
    } else {
      scale / shape * ((1 - probability)^(-shape) - 1)
    }
  }

  gpd_survival <- function(z, scale, shape) {
    if (abs(shape) < 1e-7) {
      exp(-z / scale)
    } else {
      pmax(1 + shape * z / scale, 0)^(-1 / shape)
    }
  }

  pooled_fit <- fit_pooled_gpd_mle(excess)

  # Table 2: descriptive statistics for exceedance excesses.
  descriptive_table <- data.frame(
    statistic = c(
      "Number of exceedances", "Mean excess", "Median excess", "Standard deviation",
      "Minimum excess", "Maximum excess", "Losses greater than 5% (count)",
      "Losses greater than 5% (share)"
    ),
    value = c(
      length(excess), mean(excess), stats::median(excess), stats::sd(excess),
      min(excess), max(excess), sum(excess > 5 - threshold),
      mean(excess > 5 - threshold)
    ),
    stringsAsFactors = FALSE
  )
  write_result(descriptive_table, "report_table2_descriptive.csv")

  # Table 3: pooled GPD maximum-likelihood estimates.
  pooled_table <- data.frame(
    n = length(excess),
    xi = pooled_fit$xi,
    sigma = pooled_fit$sigma,
    implied_mean_excess = if (pooled_fit$xi < 1) {
      pooled_fit$sigma / (1 - pooled_fit$xi)
    } else {
      NA_real_
    },
    loglik = pooled_fit$loglik
  )
  write_result(pooled_table, "report_table3_pooled_gpd.csv")

  # Table 4: sensitivity of the pooled GPD estimates to the threshold.
  returns <- merged$nifty_return[is.finite(merged$nifty_return)]
  threshold_grid <- c(1.5, 1.75, 2, 2.25, 2.5, 2.75, 3)
  threshold_table <- do.call(rbind, lapply(threshold_grid, function(value) {
    selected_returns <- returns[returns < -value]
    selected_excess <- -selected_returns - value

    if (length(selected_excess) < 2L) {
      return(data.frame(
        threshold = value,
        exceedances = length(selected_excess),
        mean_excess = if (length(selected_excess)) mean(selected_excess) else NA_real_,
        xi_mle = NA_real_,
        sigma_mle = NA_real_
      ))
    }

    fit <- tryCatch(
      fit_pooled_gpd_mle(selected_excess),
      error = function(e) NULL
    )

    data.frame(
      threshold = value,
      exceedances = length(selected_excess),
      mean_excess = mean(selected_excess),
      xi_mle = if (is.null(fit)) NA_real_ else fit$xi,
      sigma_mle = if (is.null(fit)) NA_real_ else fit$sigma
    )
  }))
  write_result(threshold_table, "report_table4_threshold_sensitivity.csv")

  # Full-sample MAP fits are used for coefficient and information-criterion
  # summaries. Rolling metrics are joined for side-by-side comparison.
  if (nrow(complete_exceedances) < 2L) {
    stop("There are too few complete exceedances for the full-sample MAP fits.")
  }

  constant_fit <- fit_constant_gpd(complete_exceedances$excess, prior_spec)
  constant_loglik <- sum(gpd_log_density(
    complete_exceedances$excess,
    rep(constant_fit$log_sigma, nrow(complete_exceedances)),
    constant_fit$xi
  ))

  information_criteria <- list(data.frame(
    model_set = "constant_gpd",
    prior = "baseline",
    k = 2L,
    loglik = constant_loglik
  ))
  coefficient_rows <- list()

  for (set_name in names(model_sets)) {
    selected_predictors <- model_sets[[set_name]]

    for (prior_name in c("cauchy", "lasso", "ridge", "g-prior")) {
      fitted <- fit_gpd_map(
        complete_exceedances$excess,
        complete_exceedances[, selected_predictors, drop = FALSE],
        prior = prior_name,
        prior_spec = prior_spec
      )

      predicted_scale <- predict(
        fitted, complete_exceedances[, selected_predictors, drop = FALSE]
      )
      log_likelihood <- sum(gpd_log_density(
        complete_exceedances$excess,
        log(predicted_scale),
        fitted$xi
      ))

      information_criteria[[length(information_criteria) + 1L]] <- data.frame(
        model_set = set_name,
        prior = prior_name,
        k = length(fitted$beta) + 1L,
        loglik = log_likelihood
      )

      coefficient_rows[[length(coefficient_rows) + 1L]] <- data.frame(
        model_set = set_name,
        prior = prior_name,
        term = c("(Intercept)", selected_predictors, "xi"),
        estimate = c(fitted$beta, fitted$xi),
        stringsAsFactors = FALSE
      )
    }
  }

  information_criteria <- do.call(rbind, information_criteria)
  sample_size <- nrow(complete_exceedances)
  information_criteria$AIC <-
    2 * information_criteria$k - 2 * information_criteria$loglik
  information_criteria$BIC <-
    information_criteria$k * log(sample_size) - 2 * information_criteria$loglik

  information_criteria <- merge(
    information_criteria,
    comparison_metrics[, c(
      "model_set", "prior", "RMSE", "MAE", "mean_conditional_log_score"
    )],
    by = c("model_set", "prior"),
    all.x = TRUE,
    sort = FALSE
  )
  write_result(
    information_criteria,
    "report_table8_aic_bic_and_rolling_metrics.csv"
  )

  coefficients <- do.call(rbind, coefficient_rows)
  coefficients$exp_estimate <- ifelse(
    coefficients$term == "xi", NA_real_, exp(coefficients$estimate)
  )
  write_result(coefficients, "report_table6_full_sample_coefficients.csv")

  # Preferred specification from the report: Model D with a Cauchy prior.
  preferred_set <- "D"
  preferred_prior <- "cauchy"
  preferred_predictors <- model_sets[[preferred_set]]
  preferred_fit <- fit_gpd_map(
    complete_exceedances$excess,
    complete_exceedances[, preferred_predictors, drop = FALSE],
    prior = preferred_prior,
    prior_spec = prior_spec
  )

  standardise_with_fit <- function(values, column_index) {
    (values - preferred_fit$means[column_index]) /
      preferred_fit$sds[column_index]
  }

  first_predictor_standardised <- standardise_with_fit(
    complete_exceedances[[preferred_predictors[1]]], 1
  )
  second_predictor_standardised <- standardise_with_fit(
    complete_exceedances[[preferred_predictors[2]]], 2
  )

  conditional_probability <- function(first_value, second_value) {
    fitted_scale <- exp(
      preferred_fit$beta[1] +
        preferred_fit$beta[2] * first_value +
        preferred_fit$beta[3] * second_value
    )
    gpd_survival(5 - threshold, fitted_scale, preferred_fit$xi)
  }

  percentiles <- stats::quantile(
    first_predictor_standardised,
    probs = c(0.1, 0.5, 0.9),
    names = FALSE
  )
  raw_percentiles <- stats::quantile(
    complete_exceedances[[preferred_predictors[1]]],
    probs = c(0.1, 0.5, 0.9),
    names = FALSE
  )
  second_predictor_median <- stats::median(second_predictor_standardised)

  conditional_probability_table <- data.frame(
    model = paste(preferred_set, preferred_prior),
    percentile = c(10, 50, 90),
    first_covariate_raw = raw_percentiles,
    first_covariate_standardised = percentiles,
    second_covariate_standardised_held_at = second_predictor_median,
    probability_loss_above_5_given_above_2 = vapply(
      percentiles,
      conditional_probability,
      numeric(1),
      second_value = second_predictor_median
    ),
    stringsAsFactors = FALSE
  )
  write_result(
    conditional_probability_table,
    "report_conditional_probabilities_preferred.csv"
  )

  # Figure 2: daily returns and threshold exceedances.
  open_figure("fig02_returns_exceedances.png", 11, 4.8)
  graphics::par(mar = c(3, 4.5, 2.5, 1))
  graphics::plot(
    merged$Date, merged$nifty_return,
    type = "l", col = "#7A8C9C", lwd = 0.5,
    xlab = "", ylab = "Daily log return (%)",
    main = sprintf(
      "Nifty 50 daily log returns; %d losses below -%g%% marked",
      nrow(exceedances), threshold
    )
  )
  graphics::abline(h = -threshold, col = "#C0392B", lty = 2, lwd = 1.5)
  graphics::points(
    exceedances$Date, exceedances$nifty_return,
    pch = 19, cex = 0.45, col = "#C0392B"
  )
  grDevices::dev.off()

  # Figure 3: excess-loss histogram and pooled fitted GPD density.
  open_figure("fig03_excess_histogram.png", 9, 5)
  graphics::par(mar = c(4.5, 4.5, 2.5, 1))
  graphics::hist(
    excess,
    breaks = seq(0, ceiling(max(excess)), by = 0.5),
    freq = FALSE,
    col = "#9EC1DE",
    border = "white",
    main = "Excess losses above the threshold",
    xlab = "Excess loss (percentage points)",
    ylab = "Density"
  )
  excess_grid <- seq(0, max(excess), length.out = 400)
  graphics::lines(
    excess_grid,
    exp(gpd_log_density(
      pmax(excess_grid, 1e-9),
      rep(log(pooled_fit$sigma), length(excess_grid)),
      pooled_fit$xi
    )),
    col = "#C0392B", lwd = 2
  )
  graphics::legend(
    "topright", bty = "n", lwd = c(10, 2),
    col = c("#9EC1DE", "#C0392B"),
    legend = c(
      sprintf("Observed excesses (n = %d)", length(excess)),
      sprintf("Pooled GPD MLE (xi = %.3f, sigma = %.3f)",
              pooled_fit$xi, pooled_fit$sigma)
    )
  )
  grDevices::dev.off()

  # Figure 4: empirical and fitted survival functions.
  open_figure("fig04_survival.png", 9, 5)
  graphics::par(mar = c(4.5, 4.5, 2.5, 1))
  sorted_excess <- sort(excess)
  empirical_survival <- 1 - (seq_along(sorted_excess) - 0.5) / length(sorted_excess)
  graphics::plot(
    sorted_excess, empirical_survival,
    log = "y", pch = 19, cex = 0.5, col = "#1F4E79",
    xlab = "Excess loss (percentage points)",
    ylab = "P(Z > z), logarithmic scale",
    main = "Empirical and pooled GPD survival functions",
    ylim = c(1e-3, 1)
  )
  graphics::lines(
    excess_grid,
    gpd_survival(excess_grid, pooled_fit$sigma, pooled_fit$xi),
    col = "#C0392B", lwd = 2
  )
  graphics::legend(
    "topright", bty = "n", pch = c(19, NA), lty = c(NA, 1),
    lwd = c(NA, 2), col = c("#1F4E79", "#C0392B"),
    legend = c(
      "Empirical plotting positions",
      "Pooled GPD maximum-likelihood fit"
    )
  )
  grDevices::dev.off()

  # Figure 5: GPD quantile-quantile plot.
  open_figure("fig05_qq.png", 6.5, 6.2)
  graphics::par(mar = c(4.5, 4.5, 2.5, 1))
  plotting_probabilities <- (seq_along(sorted_excess) - 0.5) / length(sorted_excess)
  theoretical_quantiles <- vapply(
    plotting_probabilities,
    gpd_quantile,
    numeric(1),
    scale = pooled_fit$sigma,
    shape = pooled_fit$xi
  )
  graphics::plot(
    theoretical_quantiles, sorted_excess,
    pch = 19, cex = 0.6, col = "#1F4E79",
    xlab = "Theoretical pooled GPD quantiles",
    ylab = "Empirical excess-loss quantiles",
    main = "GPD Q-Q plot"
  )
  graphics::abline(0, 1, lty = 2, col = "#C0392B", lwd = 1.5)
  grDevices::dev.off()

  # Figure 6: sensitivity of exceedance counts and GPD estimates to threshold.
  open_figure("fig06_threshold_sensitivity.png", 12, 4)
  graphics::par(mfrow = c(1, 3), mar = c(4.5, 4.5, 2.5, 1))
  graphics::plot(
    threshold_table$threshold, threshold_table$exceedances,
    type = "b", pch = 19, xlab = "Threshold (%)", ylab = "Exceedances",
    main = "Exceedance count"
  )
  graphics::abline(v = threshold, lty = 2, col = "grey40")
  graphics::plot(
    threshold_table$threshold, threshold_table$mean_excess,
    type = "b", pch = 19, xlab = "Threshold (%)", ylab = "Mean excess (pp)",
    main = "Mean excess"
  )
  graphics::abline(v = threshold, lty = 2, col = "grey40")
  graphics::plot(
    threshold_table$threshold, threshold_table$xi_mle,
    type = "b", pch = 19, xlab = "Threshold (%)", ylab = expression(hat(xi)),
    main = "Pooled shape estimate"
  )
  graphics::abline(v = threshold, lty = 2, col = "grey40")
  grDevices::dev.off()

  # Figure 7: loss severity against lagged realised volatility.
  open_figure("fig07_loss_vs_volatility.png", 12, 4.3)
  graphics::par(mfrow = c(1, 3), mar = c(4.5, 4.5, 2.5, 1))
  loss <- -exceedances$nifty_return
  volatility_panels <- list(
    c("nifty_vol_lag1", "(a) Nifty 50"),
    c("sp500_vol_lag1", "(b) S&P 500"),
    c("gold_vol_lag1", "(c) Gold")
  )
  for (panel in volatility_panels) {
    graphics::plot(
      exceedances[[panel[1]]] * sqrt(250), loss,
      pch = 19, cex = 0.6, col = "#2F6690",
      xlab = "Prior-date empirical volatility (annualised, %) ",
      ylab = "Loss (-log return, %)",
      main = panel[2]
    )
    graphics::abline(h = c(2, 5), lty = 2, col = c("#D97706", "#C0392B"))
  }
  grDevices::dev.off()

  # Figure 8: correlation heatmap on the complete five-predictor sample.
  correlation_matrix <- stats::cor(
    complete_exceedances[, full_predictors, drop = FALSE]
  )
  predictor_short_labels <- c(
    "Nifty vol", "S&P 500 vol", "Gold vol", "India VIX", "USD/INR vol"
  )
  open_figure("fig08_correlation_heatmap.png", 7, 6.3)
  graphics::par(mar = c(6, 7, 3, 1))
  palette <- grDevices::colorRampPalette(
    c("#2166AC", "white", "#B2182B")
  )(201)
  graphics::image(
    1:5, 1:5, correlation_matrix[, 5:1],
    col = palette, zlim = c(-1, 1), axes = FALSE,
    xlab = "", ylab = "",
    main = sprintf(
      "Pearson correlations (%d complete exceedances)",
      nrow(complete_exceedances)
    )
  )
  graphics::axis(1, 1:5, predictor_short_labels, las = 2)
  graphics::axis(2, 1:5, rev(predictor_short_labels), las = 1)
  for (i in 1:5) {
    for (j in 1:5) {
      graphics::text(i, 6 - j, sprintf("%.2f", correlation_matrix[i, j]), cex = 0.95)
    }
  }
  grDevices::dev.off()

  # Figure 9: number of exceedances in each rolling training window.
  open_figure("fig09_training_exceedances.png", 10, 4.5)
  graphics::par(mar = c(3, 4.5, 2.5, 1))
  graphics::plot(
    comparison_windows$forecast_date,
    comparison_windows$n_training_exceedances,
    type = "b", pch = 19, col = "#1F4E79",
    xlab = "", ylab = "Training exceedances in window",
    main = "Training exceedances available at each forecast date"
  )
  graphics::abline(
    h = config$minimum_training_exceedances,
    lty = 2, col = "#C0392B"
  )
  grDevices::dev.off()

  get_model_predictions <- function(model_set, prior) {
    data <- comparison_predictions[
      comparison_predictions$model_set == model_set &
        comparison_predictions$prior == prior,
      , drop = FALSE
    ]
    data[order(data$Date), , drop = FALSE]
  }

  selected_models <- list(
    D = get_model_predictions("D", "cauchy"),
    C = get_model_predictions("C", "g-prior"),
    full = get_model_predictions("full_5", "g-prior")
  )
  model_colours <- c(
    constant = "#7F8C8D",
    C = "#2563EB",
    D = "#D97706",
    full = "#2E8B57"
  )
  model_labels <- c(
    D = "Model D (Cauchy)",
    C = "Model C (g-prior)",
    full = "Full five-predictor model (g-prior)"
  )

  # Figures 10 and 11: rolling shape and scale estimates.
  for (quantity in c("shape", "predicted_scale")) {
    filename <- if (quantity == "shape") {
      "fig10_rolling_shape.png"
    } else {
      "fig11_rolling_scale.png"
    }
    y_range <- range(unlist(lapply(selected_models, function(model) model[[quantity]])), finite = TRUE)

    open_figure(filename, 10, 4.5)
    graphics::par(mar = c(3, 4.5, 2.5, 1))
    graphics::plot(
      selected_models$D$Date,
      selected_models$D[[quantity]],
      type = "n", ylim = y_range,
      xlab = "",
      ylab = if (quantity == "shape") {
        expression(hat(xi))
      } else {
        expression(hat(sigma) ~ "(percentage points)")
      },
      main = if (quantity == "shape") {
        "Rolling GPD shape estimates"
      } else {
        "Rolling predicted GPD scale"
      }
    )
    if (quantity == "shape") graphics::abline(h = 0, lty = 3)

    for (model_name in names(selected_models)) {
      graphics::points(
        selected_models[[model_name]]$Date,
        selected_models[[model_name]][[quantity]],
        pch = 19, cex = 0.8,
        col = model_colours[model_name]
      )
    }
    graphics::legend(
      "topleft", bty = "n", pch = 19,
      col = model_colours[names(selected_models)],
      legend = model_labels[names(selected_models)]
    )
    grDevices::dev.off()
  }

  # Figure 12: preferred model fitted values and conditional tail probabilities.
  preferred_fitted_mean <- predict(
    preferred_fit,
    complete_exceedances[, preferred_predictors, drop = FALSE]
  ) / (1 - preferred_fit$xi)
  first_predictor_grid <- seq(
    min(first_predictor_standardised),
    max(first_predictor_standardised),
    length.out = 200
  )
  probability_grid <- vapply(
    first_predictor_grid,
    conditional_probability,
    numeric(1),
    second_value = second_predictor_median
  )

  open_figure("fig12_preferred_model_fit.png", 11, 4.8)
  graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1))
  graphics::plot(
    preferred_fitted_mean, complete_exceedances$excess,
    pch = 19, cex = 0.6, col = "#D97706",
    xlab = "Fitted conditional mean excess (pp)",
    ylab = "Observed excess (pp)",
    main = "(a) Observed and fitted mean excess"
  )
  graphics::abline(0, 1, lty = 2)

  graphics::plot(
    first_predictor_grid, probability_grid,
    type = "l", lwd = 2, col = "#D97706",
    ylim = c(0, max(probability_grid, na.rm = TRUE)),
    xlab = paste("Standardised prior-date", "India VIX"),
    ylab = "P(L > 5% | L > 2%, x)",
    main = "(b) Conditional probability of a loss above 5%"
  )
  graphics::abline(v = percentiles, lty = 3, col = "grey40")
  graphics::text(
    percentiles, rep(0, length(percentiles)),
    labels = c("P10", "P50", "P90"), pos = 4, cex = 0.8
  )
  grDevices::dev.off()

  # Figure 13: observed excesses and selected forecasts on the common dates.
  constant_predictions <- get_model_predictions("constant_gpd", "baseline")
  open_figure("fig13_obs_vs_pred.png", 11, 4.8)
  graphics::par(mar = c(3, 4.5, 2.5, 1))
  graphics::plot(
    constant_predictions$Date,
    constant_predictions$observed_excess,
    pch = 19, col = "black",
    ylim = c(0, max(constant_predictions$observed_excess, na.rm = TRUE) * 1.05),
    xlab = "", ylab = "Excess loss (pp)",
    main = "Observed excess and predicted conditional mean (common forecast dates)"
  )
  graphics::points(
    constant_predictions$Date,
    constant_predictions$predicted_mean_excess,
    pch = 15, col = model_colours["constant"]
  )
  graphics::points(
    selected_models$C$Date,
    selected_models$C$predicted_mean_excess,
    pch = 17, col = model_colours["C"]
  )
  graphics::points(
    selected_models$D$Date,
    selected_models$D$predicted_mean_excess,
    pch = 18, cex = 1.2, col = model_colours["D"]
  )
  graphics::legend(
    "topleft", bty = "n", pch = c(19, 15, 17, 18),
    col = c("black", model_colours[c("constant", "C", "D")]),
    legend = c(
      "Observed", "Constant GPD", "Model C (g-prior)", "Model D (Cauchy)"
    )
  )
  grDevices::dev.off()

  # Figure 14: RMSE, MAE, and conditional log scores for all prior families.
  open_figure("fig14_metrics_by_set.png", 12, 4.3)
  graphics::par(mfrow = c(1, 3), mar = c(4.5, 4.5, 2.5, 1))
  set_order <- c("full_5", "A", "B", "C", "D")
  prior_order <- c("cauchy", "lasso", "ridge", "g-prior")
  prior_colours <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A")
  offsets <- c(-0.15, -0.05, 0.05, 0.15)

  for (metric_name in c("RMSE", "MAE", "mean_conditional_log_score")) {
    subset_metrics <- comparison_metrics[
      comparison_metrics$model_set %in% set_order,
      , drop = FALSE
    ]
    baseline_value <- comparison_metrics[
      comparison_metrics$model_set == "constant_gpd",
      metric_name
    ]

    metric_values <- subset_metrics[[metric_name]]
    metric_range <- range(c(metric_values, baseline_value), finite = TRUE)
    graphics::plot(
      NA,
      xlim = c(0.6, 5.4),
      ylim = metric_range,
      xaxt = "n",
      xlab = "Predictor set",
      ylab = metric_name,
      main = if (metric_name == "mean_conditional_log_score") {
        "Mean conditional log score (higher is better)"
      } else {
        paste(metric_name, "(lower is better)")
      }
    )
    graphics::axis(1, 1:5, c("Full-5", "A", "B", "C", "D"))
    graphics::abline(h = baseline_value, lty = 3)

    for (j in seq_along(prior_order)) {
      prior_data <- subset_metrics[subset_metrics$prior == prior_order[j], ]
      x_position <- match(prior_data$model_set, set_order) + offsets[j]
      graphics::points(
        x_position, prior_data[[metric_name]],
        pch = 19, col = prior_colours[j]
      )
    }

    if (metric_name == "RMSE") {
      graphics::legend(
        "topright", bty = "n", pch = c(rep(19, 4), NA),
        lty = c(rep(NA, 4), 3),
        col = c(prior_colours, "black"),
        legend = c(prior_order, "Constant GPD baseline"),
        cex = 0.85
      )
    }
  }
  grDevices::dev.off()

  # Figure 15: probability integral transform (PIT) diagnostics.
  calculate_pit <- function(data) {
    1 - mapply(
      gpd_survival,
      data$observed_excess,
      data$predicted_scale,
      data$shape
    )
  }

  pit_values <- list(
    constant = calculate_pit(constant_predictions),
    C = calculate_pit(selected_models$C),
    D = calculate_pit(selected_models$D)
  )
  ks_p_values <- vapply(pit_values, function(values) {
    suppressWarnings(stats::ks.test(values, "punif")$p.value)
  }, numeric(1))

  pit_summary <- data.frame(
    model = names(pit_values),
    ks_p_value = ks_p_values,
    median_pit = vapply(pit_values, stats::median, numeric(1)),
    stringsAsFactors = FALSE
  )
  write_result(pit_summary, "report_pit_summary.csv")

  open_figure("fig15_pit.png", 11, 4.8)
  graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1))
  graphics::plot(
    c(0, 1), c(0, 1), type = "l", lty = 2, col = "grey40",
    xlab = "PIT value", ylab = "Empirical CDF",
    main = "PIT empirical CDF"
  )
  pit_colours <- c(
    constant = model_colours["constant"],
    C = model_colours["C"],
    D = model_colours["D"]
  )
  for (model_name in names(pit_values)) {
    graphics::lines(
      stats::ecdf(pit_values[[model_name]]),
      col = pit_colours[model_name],
      do.points = FALSE,
      verticals = TRUE,
      lwd = 1.6
    )
  }
  graphics::legend(
    "topleft", bty = "n", lwd = 1.6,
    col = pit_colours[names(pit_values)],
    legend = sprintf(
      "%s (KS p = %.2f)",
      c("Constant GPD", "Model C (g-prior)", "Model D (Cauchy)"),
      ks_p_values
    )
  )

  graphics::plot(
    c(0, 1), c(0, 1), type = "l", lty = 2, col = "grey40",
    xlab = "Uniform quantile", ylab = "Ordered PIT",
    main = "PIT Q-Q plot"
  )
  for (model_name in names(pit_values)) {
    values <- pit_values[[model_name]]
    graphics::points(
      (seq_along(values) - 0.5) / length(values),
      sort(values),
      pch = 19, cex = 0.6,
      col = pit_colours[model_name]
    )
  }
  grDevices::dev.off()

  message("Report tables saved to: ", results_dir)
  message("Report figures saved to: ", figure_dir)

  invisible(list(
    descriptive = descriptive_table,
    pooled_gpd = pooled_table,
    threshold_sensitivity = threshold_table,
    information_criteria = information_criteria,
    coefficients = coefficients,
    conditional_probabilities = conditional_probability_table,
    pit_summary = pit_summary
  ))
}
