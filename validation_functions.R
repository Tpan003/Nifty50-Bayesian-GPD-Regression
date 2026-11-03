# Checks for data preparation, likelihood calculations, rolling forecasts,
# and the consistency of the exported model metrics.

validate_analysis <- function(merged_data, predictions, diagnostics, metrics,
                              input_path, config) {
  checks <- list()

  add_check <- function(name, passed, detail) {
    checks[[length(checks) + 1L]] <<- data.frame(
      check = name,
      passed = isTRUE(passed),
      detail = as.character(detail),
      stringsAsFactors = FALSE
    )
  }

  predictors <- c(
    "nifty_vol_lag1", "sp500_vol_lag1", "gold_vol_lag1",
    "india_vix_lag1", "usdinr_vol_lag1"
  )

  # Dates and exceedance construction.
  dates_ok <- !anyNA(merged_data$Date) &&
    !anyDuplicated(merged_data$Date) &&
    identical(order(merged_data$Date), seq_len(nrow(merged_data)))
  add_check(
    "Nifty dates are unique and sorted",
    dates_ok,
    sprintf("%d forecast dates", nrow(merged_data))
  )

  expected_event <- as.integer(
    is.finite(merged_data$nifty_return) &
      merged_data$nifty_return < -config$threshold
  )
  response_ok <- identical(merged_data$exceedance, expected_event) &&
    all(merged_data$excess[expected_event == 1L] > 0) &&
    all(is.na(merged_data$excess[expected_event == 0L]))
  add_check(
    "Threshold indicator and excess response are correct",
    response_ok,
    sprintf(
      "Threshold %.2f percentage points; %d exceedances",
      config$threshold, sum(expected_event)
    )
  )

  # A small holiday/weekend example checks that same-day observations are not used.
  toy_target <- as.Date(c("2024-01-08", "2024-01-09", "2024-01-10"))
  toy_source <- as.Date(c("2024-01-05", "2024-01-08", "2024-01-12"))
  toy_result <- latest_prior(toy_target, toy_source, c(5, 8, 12))
  add_check(
    "Lag function uses strictly prior observations",
    identical(unname(toy_result), c(5, 8, 8)),
    "An observation on the target date is excluded; non-trading days use the previous available observation."
  )

  # Compare the lagged columns with a separate direct lookup on sampled dates.
  source_sheets <- read_market_sheets(input_path)
  get_date_column <- function(sheet) {
    if ("DATE" %in% names(sheet)) "DATE" else "Date"
  }
  source_map <- list(
    nifty_vol_lag1 = list(
      sheet = source_sheets$nifty,
      value = "Daily empirical volatility (%)"
    ),
    sp500_vol_lag1 = list(
      sheet = source_sheets$sp500,
      value = "Daily empirical volatility (%)"
    ),
    gold_vol_lag1 = list(
      sheet = source_sheets$gold,
      value = "Daily empirical volatility (%)"
    ),
    india_vix_lag1 = list(
      sheet = source_sheets$vix,
      value = "India VIX level (points)"
    ),
    usdinr_vol_lag1 = list(
      sheet = source_sheets$usdinr,
      value = "Daily empirical volatility (%)"
    )
  )

  sample_rows <- unique(as.integer(round(seq(
    1, nrow(merged_data), length.out = min(100L, nrow(merged_data))
  ))))

  for (column_name in names(source_map)) {
    source_info <- source_map[[column_name]]
    source_sheet <- source_info$sheet
    source_dates <- as.Date(source_sheet[[get_date_column(source_sheet)]])
    source_values <- as.numeric(source_sheet[[source_info$value]])

    matches <- vapply(sample_rows, function(row_index) {
      available <- which(
        !is.na(source_dates) &
          source_dates < merged_data$Date[row_index] &
          is.finite(source_values)
      )

      if (length(available) == 0L) {
        return(is.na(merged_data[[column_name]][row_index]))
      }

      latest_row <- available[which.max(source_dates[available])]
      is.finite(merged_data[[column_name]][row_index]) &&
        abs(merged_data[[column_name]][row_index] -
              source_values[latest_row]) < 1e-10
    }, logical(1))

    add_check(
      paste("Actual prior-date alignment:", column_name),
      all(matches),
      sprintf("Independent lookup on %d sampled Nifty dates", length(sample_rows))
    )
  }

  # Likelihood checks: exponential limit, numerical integration, and support.
  test_z <- c(0.1, 0.7, 1.8)
  test_sigma <- 1.4
  exponential_limit <- gpd_log_density(
    test_z, rep(log(test_sigma), length(test_z)), 0
  )
  exact_exponential <- stats::dexp(
    test_z, rate = 1 / test_sigma, log = TRUE
  )
  add_check(
    "GPD exponential limit",
    isTRUE(all.equal(exponential_limit, exact_exponential, tolerance = 1e-12)),
    "At xi = 0, the GPD density agrees with the exponential density."
  )

  density_integral <- tryCatch(
    stats::integrate(
      function(z) exp(gpd_log_density(
        z, rep(log(test_sigma), length(z)), 0.2
      )),
      lower = 0, upper = Inf, rel.tol = 1e-8
    )$value,
    error = function(e) NA_real_
  )
  add_check(
    "Positive-shape GPD density integrates to one",
    is.finite(density_integral) && abs(density_integral - 1) < 1e-7,
    sprintf("Numerical integral: %.10f", density_integral)
  )

  support_ok <-
    is.infinite(gpd_log_density(2.1, log(1), -0.5)) &&
    gpd_log_density(2.1, log(1), -0.5) < 0 &&
    is.finite(gpd_log_density(1, log(1), -0.5))
  add_check(
    "GPD finite support is enforced",
    support_ok,
    "For xi = -0.5 and sigma = 1, z = 1 is valid and z = 2.1 is outside support."
  )

  # Compare analytic gradients with central finite differences.
  gradient_x <- cbind(
    Intercept = rep(1, 5),
    x1 = c(-1, -0.5, 0, 0.5, 1),
    x2 = c(0.5, -1, 0, 1, -0.5),
    x3 = c(1, 0.2, -0.7, 0.3, -0.8)
  )
  gradient_z <- c(0.2, 0.35, 0.7, 1.1, 0.5)
  gradient_parameters <- c(log(0.8), 0.2, -0.15, 0.1, 0.18)
  gradient_errors <- vapply(
    c("cauchy", "lasso", "ridge", "g-prior"),
    function(prior) {
      analytic <- gpd_map_gradient(
        gradient_parameters, gradient_z, gradient_x, prior
      )
      numeric <- vapply(seq_along(gradient_parameters), function(j) {
        step <- 1e-6
        plus <- minus <- gradient_parameters
        plus[j] <- plus[j] + step
        minus[j] <- minus[j] - step
        (
          gpd_map_objective(plus, gradient_z, gradient_x, prior) -
            gpd_map_objective(minus, gradient_z, gradient_x, prior)
        ) / (2 * step)
      }, numeric(1))
      max(abs(analytic - numeric))
    },
    numeric(1)
  )
  add_check(
    "Analytic gradients agree with finite differences",
    all(is.finite(gradient_errors)) && max(gradient_errors) < 1e-4,
    sprintf("Maximum absolute difference: %.3g", max(gradient_errors))
  )

  # Confirm that each prior is centred at zero and that the g-prior uses X'X/(2g).
  prior_design <- cbind(
    Intercept = 1,
    x1 = c(-1, 0, 1, 1, 0, -1),
    x2 = c(1, -1, 0, 0, 1, -1)
  )
  zero_penalties <- vapply(
    c("cauchy", "lasso", "ridge", "g-prior"),
    function(prior) {
      prior_penalty(c(0, 0, 0), prior, prior_design, config$prior_specification)
    },
    numeric(1)
  )
  test_beta <- c(0, 0.2, -0.1)
  g_scale <- if (identical(config$prior_specification$g_prior_g, "n")) {
    nrow(prior_design)
  } else {
    as.numeric(config$prior_specification$g_prior_g)
  }
  observed_g_penalty <- prior_penalty(
    test_beta, "g-prior", prior_design, config$prior_specification
  )
  expected_g_penalty <- as.numeric(
    crossprod(
      test_beta[-1L],
      crossprod(prior_design[, -1L, drop = FALSE]) %*% test_beta[-1L]
    )
  ) / (2 * g_scale)
  add_check(
    "Prior penalties are correctly centred and specified",
    all(abs(zero_penalties) < 1e-12) &&
      isTRUE(all.equal(observed_g_penalty, expected_g_penalty, tolerance = 1e-12)),
    "All slope penalties are zero at zero; the g-prior uses X'X divided by 2g."
  )

  # Forecast presence and forecast-row consistency.
  has_predictions <- nrow(predictions) > 0L
  add_check(
    "Rolling predictions were produced",
    has_predictions,
    sprintf("%d model-date rows", nrow(predictions))
  )

  if (has_predictions) {
    row_index <- match(predictions$Date, merged_data$Date)
    event_rows_ok <- !anyNA(row_index) &&
      all(merged_data$exceedance[row_index] == 1L) &&
      isTRUE(all.equal(
        predictions$observed_excess,
        merged_data$excess[row_index],
        tolerance = 1e-12
      )) &&
      all(complete.cases(merged_data[row_index, predictors, drop = FALSE]))

    add_check(
      "Forecast responses match eligible exceedance rows",
      event_rows_ok,
      "Every scored response matches its source row and has complete lagged predictors."
    )

    scales_ok <- all(
      is.finite(predictions$predicted_scale) &
        predictions$predicted_scale > 0
    )
    add_check(
      "Predicted GPD scales are positive and finite",
      scales_ok,
      sprintf("%d forecasts checked", nrow(predictions))
    )

    score_ok <- !anyNA(predictions$log_score) &&
      !any(is.nan(predictions$log_score)) &&
      !any(predictions$log_score == Inf)
    add_check(
      "Held-out conditional log scores are valid",
      score_ok,
      sprintf(
        "%d predictions received -Inf because observations were outside fitted support",
        sum(predictions$log_score == -Inf)
      )
    )

    regression_models <- c("cauchy", "lasso", "ridge", "g-prior")
    regression_predictions <- predictions[
      predictions$model %in% regression_models,
      , drop = FALSE
    ]
    prediction_groups <- split(regression_predictions, regression_predictions$model)
    distinct_scales <- all(regression_models %in% names(prediction_groups))
    if (distinct_scales) {
      common_dates <- Reduce(
        intersect,
        lapply(prediction_groups[regression_models], function(data) data$Date)
      )
      scale_matrix <- do.call(cbind, lapply(regression_models, function(model) {
        data <- prediction_groups[[model]]
        data$predicted_scale[match(common_dates, data$Date)]
      }))
      distinct_scales <- length(common_dates) > 0L &&
        any(apply(scale_matrix, 1L, function(values) {
          max(values) - min(values) > 1e-8
        }))
    }
    add_check(
      "Prior-specific models produce distinct predicted scales",
      distinct_scales,
      "At least one common forecast date has different scales across the four prior models."
    )
  }

  if (nrow(diagnostics) > 0L) {
    unique_fit_keys <- !anyDuplicated(
      paste(diagnostics$Date, diagnostics$model, sep = "|")
    )
    add_check(
      "One diagnostic row exists per date and model",
      unique_fit_keys,
      sprintf("%d attempted fits", nrow(diagnostics))
    )

    successful <- !is.na(diagnostics$convergence) & diagnostics$convergence == 0L
    successful_fit_details_ok <-
      all(is.finite(diagnostics$objective[successful])) &&
      all(nzchar(diagnostics$message[successful]))
    add_check(
      "Successful fits have finite objectives and status messages",
      successful_fit_details_ok,
      sprintf(
        "%d successful; %d failed or skipped",
        sum(successful), sum(!successful)
      )
    )

    regression_models <- c("cauchy", "lasso", "ridge", "g-prior")
    nonzero_slope_evidence <- vapply(regression_models, function(model) {
      rows <- diagnostics$model == model & successful
      any(
        rows & is.finite(diagnostics$max_abs_slope) &
          diagnostics$max_abs_slope > 1e-8
      )
    }, logical(1))
    add_check(
      "Prior-regularised regressions estimate non-zero slopes",
      all(nonzero_slope_evidence),
      paste(
        names(nonzero_slope_evidence), nonzero_slope_evidence,
        collapse = "; "
      )
    )

    window_checks <- vapply(seq_len(nrow(diagnostics)), function(i) {
      forecast_date <- diagnostics$Date[i]
      earlier_dates <- sort(unique(
        merged_data$Date[merged_data$Date < forecast_date]
      ))
      window_dates <- tail(earlier_dates, config$rolling_observations)
      count <- sum(
        merged_data$Date %in% window_dates &
          merged_data$exceedance == 1L &
          is.finite(merged_data$excess) &
          merged_data$excess > 0 &
          complete.cases(merged_data[, predictors, drop = FALSE])
      )
      all(window_dates < forecast_date) && count == diagnostics$n_train[i]
    }, logical(1))

    add_check(
      "Rolling windows use strictly prior dates and correct event counts",
      all(window_checks),
      sprintf("%d diagnostic rows independently checked", length(window_checks))
    )

    # Refit the first ridge fold and independently verify that prediction uses
    # the means and standard deviations learned from that training fold only.
    first_forecast_date <- min(diagnostics$Date)
    earlier_dates <- sort(unique(
      merged_data$Date[merged_data$Date < first_forecast_date]
    ))
    training_dates <- tail(earlier_dates, config$rolling_observations)
    training_window <- merged_data[
      merged_data$Date %in% training_dates,
      , drop = FALSE
    ]
    training_events <- training_window[
      training_window$exceedance == 1L &
        is.finite(training_window$excess) &
        training_window$excess > 0 &
        complete.cases(training_window[, predictors, drop = FALSE]),
      , drop = FALSE
    ]
    target <- merged_data[
      merged_data$Date == first_forecast_date &
        merged_data$exceedance == 1L &
        complete.cases(merged_data[, predictors, drop = FALSE]),
      , drop = FALSE
    ]

    scaling_ok <- FALSE
    if (nrow(training_events) >= 2L && nrow(target) > 0L) {
      probe <- tryCatch(
        fit_gpd_map(
          training_events$excess,
          training_events[, predictors, drop = FALSE],
          "ridge",
          config$prior_specification
        ),
        error = function(error) NULL
      )
      if (!is.null(probe)) {
        manual_design <- cbind(
          Intercept = 1,
          scale(
            as.matrix(target[, predictors, drop = FALSE]),
            center = probe$means,
            scale = probe$sds
          )
        )
        manual_scale <- exp(drop(manual_design %*% probe$beta))
        scaling_ok <-
          isTRUE(all.equal(
            unname(probe$means),
            unname(colMeans(training_events[, predictors, drop = FALSE])),
            tolerance = 1e-12
          )) &&
          isTRUE(all.equal(
            unname(probe$sds),
            unname(apply(training_events[, predictors, drop = FALSE], 2, stats::sd)),
            tolerance = 1e-12
          )) &&
          isTRUE(all.equal(
            unname(predict(probe, target[, predictors, drop = FALSE])),
            unname(manual_scale),
            tolerance = 1e-12
          ))
      }
    }
    add_check(
      "Training-only standardisation is reused for prediction",
      scaling_ok,
      "The first rolling ridge fold was reconstructed and checked independently."
    )
  } else {
    add_check(
      "Rolling fit diagnostics are present",
      FALSE,
      "No eligible rolling fits were recorded."
    )
  }

  # Recalculate RMSE, MAE, and mean conditional log score from prediction rows.
  if (nrow(metrics) > 0L && has_predictions) {
    metric_checks <- vapply(seq_len(nrow(metrics)), function(i) {
      model_data <- predictions[
        predictions$model == metrics$model[i], , drop = FALSE
      ]
      valid_mean <- is.finite(model_data$predicted_mean_excess)
      calculated_rmse <- if (any(valid_mean)) {
        sqrt(mean((model_data$observed_excess[valid_mean] -
                     model_data$predicted_mean_excess[valid_mean])^2))
      } else {
        NA_real_
      }
      calculated_mae <- if (any(valid_mean)) {
        mean(abs(model_data$observed_excess[valid_mean] -
                   model_data$predicted_mean_excess[valid_mean]))
      } else {
        NA_real_
      }

      metrics$forecasts[i] == nrow(model_data) &&
        isTRUE(all.equal(metrics$RMSE[i], calculated_rmse, tolerance = 1e-12)) &&
        isTRUE(all.equal(metrics$MAE[i], calculated_mae, tolerance = 1e-12)) &&
        isTRUE(all.equal(
          metrics$mean_conditional_log_score[i],
          mean(model_data$log_score),
          tolerance = 1e-12
        ))
    }, logical(1))

    add_check(
      "Metric table matches the prediction rows",
      all(metric_checks),
      sprintf("%d model summaries independently checked", length(metric_checks))
    )
  } else {
    add_check("Metric table is present", FALSE, "No model metrics were generated.")
  }

  result <- do.call(rbind, checks)
  rownames(result) <- NULL

  if (any(!result$passed)) {
    failed_checks <- result$check[!result$passed]
    stop(
      "Validation failed: ",
      paste(failed_checks, collapse = "; ")
    )
  }

  result
}
