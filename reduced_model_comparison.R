# Compare the full five-predictor model with prespecified reduced predictor
# sets. Every model is evaluated on the same forecast dates and training folds.

reduced_model_sets <- function() {
  list(
    full_5 = c(
      "nifty_vol_lag1", "sp500_vol_lag1", "gold_vol_lag1",
      "india_vix_lag1", "usdinr_vol_lag1"
    ),
    A = c(
      "nifty_vol_lag1", "sp500_vol_lag1", "gold_vol_lag1",
      "usdinr_vol_lag1"
    ),
    B = c(
      "sp500_vol_lag1", "gold_vol_lag1", "india_vix_lag1",
      "usdinr_vol_lag1"
    ),
    C = c("nifty_vol_lag1", "usdinr_vol_lag1"),
    D = c("india_vix_lag1", "usdinr_vol_lag1")
  )
}

reduced_design_diagnostics <- function(x) {
  x <- as.data.frame(x)
  if (ncol(x) < 2L) {
    return(c(max_abs_pairwise_r = NA_real_, max_vif = 1))
  }

  correlation <- stats::cor(x)
  max_abs_r <- max(abs(correlation[upper.tri(correlation)]))

  vif <- vapply(names(x), function(response) {
    other_predictors <- setdiff(names(x), response)
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

  c(max_abs_pairwise_r = max_abs_r, max_vif = max(vif))
}

run_reduced_model_comparison <- function(
    merged_data,
    config,
    baseline_predictions,
    baseline_diagnostics,
    output_dir = "results") {

  model_sets <- reduced_model_sets()
  full_predictors <- model_sets$full_5
  priors <- c("cauchy", "lasso", "ridge", "g-prior")
  required_columns <- c("Date", "exceedance", "excess", full_predictors)

  if (!all(required_columns %in% names(merged_data))) {
    stop("Merged data is missing columns required for model comparison.")
  }
  if (is.null(baseline_predictions) || nrow(baseline_predictions) == 0L ||
      is.null(baseline_diagnostics) || nrow(baseline_diagnostics) == 0L) {
    stop("Pass the original rolling predictions and diagnostics to define the common folds.")
  }

  merged_data$Date <- as.Date(merged_data$Date)
  forecast_dates <- sort(unique(as.Date(baseline_predictions$Date)))
  if (length(forecast_dates) == 0L || anyNA(forecast_dates)) {
    stop("No valid forecast dates were supplied by the original rolling pipeline.")
  }

  complete_full <- stats::complete.cases(
    merged_data[, full_predictors, drop = FALSE]
  )
  eligible_dates <- merged_data$Date[
    merged_data$exceedance == 1L &
      is.finite(merged_data$excess) &
      merged_data$excess > 0 &
      complete_full
  ]
  if (!all(forecast_dates %in% eligible_dates)) {
    stop("At least one requested forecast date is not a complete exceedance.")
  }

  prediction_rows <- list()
  fit_rows <- list()
  window_rows <- list()

  for (date_value in forecast_dates) {
    # Restore Date class because base R can drop it in a for-loop.
    forecast_date <- as.Date(date_value, origin = "1970-01-01")

    previous_dates <- sort(unique(
      merged_data$Date[merged_data$Date < forecast_date]
    ))
    window_dates <- tail(previous_dates, config$rolling_observations)

    if (length(window_dates) != config$rolling_observations ||
        max(window_dates) >= forecast_date) {
      stop("A comparison fold does not contain the configured prior-date window.")
    }

    training_window <- merged_data[
      merged_data$Date %in% window_dates, , drop = FALSE
    ]
    complete_window <- stats::complete.cases(
      training_window[, full_predictors, drop = FALSE]
    )

    training_events <- training_window[
      training_window$exceedance == 1L &
        is.finite(training_window$excess) &
        training_window$excess > 0 &
        complete_window,
      , drop = FALSE
    ]
    target <- merged_data[
      merged_data$Date == forecast_date &
        merged_data$exceedance == 1L &
        is.finite(merged_data$excess) &
        merged_data$excess > 0 &
        complete_full,
      , drop = FALSE
    ]

    if (nrow(training_events) < config$minimum_training_exceedances ||
        nrow(target) != 1L) {
      stop("A comparison fold does not meet the training-event or target-row requirement.")
    }

    window_rows[[length(window_rows) + 1L]] <- data.frame(
      forecast_date = forecast_date,
      window_first_date = min(window_dates),
      window_last_date = max(window_dates),
      n_training_exceedances = nrow(training_events),
      n_window_dates = length(window_dates),
      stringsAsFactors = FALSE
    )

    # Intercept-only benchmark, evaluated on the same fold as every regression.
    baseline_fit <- tryCatch(
      fit_constant_gpd(training_events$excess, config$prior_specification),
      error = function(e) e
    )
    if (inherits(baseline_fit, "error")) {
      stop(
        "Constant GPD fit failed on ", forecast_date, ": ",
        conditionMessage(baseline_fit)
      )
    }

    baseline_scale <- exp(baseline_fit$log_sigma)
    baseline_mean <- if (baseline_fit$xi < 1) {
      baseline_scale / (1 - baseline_fit$xi)
    } else {
      NA_real_
    }

    prediction_rows[[length(prediction_rows) + 1L]] <- data.frame(
      Date = forecast_date,
      model_set = "constant_gpd",
      prior = "baseline",
      predictor_count = 0L,
      predictors = "(intercept only)",
      observed_excess = target$excess,
      predicted_scale = baseline_scale,
      predicted_mean_excess = baseline_mean,
      shape = baseline_fit$xi,
      log_score = gpd_log_density(
        target$excess, baseline_fit$log_sigma, baseline_fit$xi
      ),
      stringsAsFactors = FALSE
    )

    fit_rows[[length(fit_rows) + 1L]] <- data.frame(
      forecast_date = forecast_date,
      model_set = "constant_gpd",
      prior = "baseline",
      n_training_exceedances = nrow(training_events),
      predictor_count = 0L,
      max_abs_pairwise_r = NA_real_,
      max_vif = NA_real_,
      convergence = baseline_fit$convergence,
      objective = baseline_fit$objective,
      iterations = baseline_fit$iterations,
      max_abs_slope = 0,
      error = "",
      stringsAsFactors = FALSE
    )

    # Fit all four prior specifications for each prespecified predictor set.
    for (set_name in names(model_sets)) {
      selected_predictors <- model_sets[[set_name]]
      training_x <- training_events[, selected_predictors, drop = FALSE]
      target_x <- target[, selected_predictors, drop = FALSE]
      design_stats <- reduced_design_diagnostics(training_x)

      for (prior_name in priors) {
        fitted <- tryCatch(
          fit_gpd_map(
            training_events$excess,
            training_x,
            prior = prior_name,
            prior_spec = config$prior_specification
          ),
          error = function(e) e
        )

        if (inherits(fitted, "error")) {
          fit_rows[[length(fit_rows) + 1L]] <- data.frame(
            forecast_date = forecast_date,
            model_set = set_name,
            prior = prior_name,
            n_training_exceedances = nrow(training_events),
            predictor_count = length(selected_predictors),
            max_abs_pairwise_r = unname(design_stats["max_abs_pairwise_r"]),
            max_vif = unname(design_stats["max_vif"]),
            convergence = NA_integer_,
            objective = NA_real_,
            iterations = NA_integer_,
            max_abs_slope = NA_real_,
            error = conditionMessage(fitted),
            stringsAsFactors = FALSE
          )
          next
        }

        predicted_scale <- predict(fitted, target_x)
        predicted_mean <- if (fitted$xi < 1) {
          predicted_scale / (1 - fitted$xi)
        } else {
          NA_real_
        }
        log_score <- gpd_log_density(
          target$excess, log(predicted_scale), fitted$xi
        )

        prediction_rows[[length(prediction_rows) + 1L]] <- data.frame(
          Date = forecast_date,
          model_set = set_name,
          prior = prior_name,
          predictor_count = length(selected_predictors),
          predictors = paste(selected_predictors, collapse = ";"),
          observed_excess = target$excess,
          predicted_scale = predicted_scale,
          predicted_mean_excess = predicted_mean,
          shape = fitted$xi,
          log_score = log_score,
          stringsAsFactors = FALSE
        )

        fit_rows[[length(fit_rows) + 1L]] <- data.frame(
          forecast_date = forecast_date,
          model_set = set_name,
          prior = prior_name,
          n_training_exceedances = nrow(training_events),
          predictor_count = length(selected_predictors),
          max_abs_pairwise_r = unname(design_stats["max_abs_pairwise_r"]),
          max_vif = unname(design_stats["max_vif"]),
          convergence = fitted$convergence,
          objective = fitted$objective,
          iterations = fitted$iterations,
          max_abs_slope = max(abs(fitted$beta[-1L])),
          error = "",
          stringsAsFactors = FALSE
        )
      }
    }
  }

  predictions <- do.call(rbind, prediction_rows)
  diagnostics <- do.call(rbind, fit_rows)
  windows <- do.call(rbind, window_rows)

  expected_sets <- c("constant_gpd", names(model_sets))
  expected_rows <- length(forecast_dates) *
    (1L + length(model_sets) * length(priors))

  groups <- split(
    predictions,
    interaction(predictions$model_set, predictions$prior, drop = TRUE)
  )
  metrics <- do.call(rbind, lapply(groups, function(data) {
    valid_mean <- is.finite(data$predicted_mean_excess)

    data.frame(
      model_set = data$model_set[1],
      prior = data$prior[1],
      predictor_count = data$predictor_count[1],
      predictors = data$predictors[1],
      forecasts = nrow(data),
      finite_mean_forecasts = sum(valid_mean),
      RMSE = if (any(valid_mean)) {
        sqrt(mean((data$observed_excess[valid_mean] -
                     data$predicted_mean_excess[valid_mean])^2))
      } else {
        NA_real_
      },
      MAE = if (any(valid_mean)) {
        mean(abs(data$observed_excess[valid_mean] -
                   data$predicted_mean_excess[valid_mean]))
      } else {
        NA_real_
      },
      mean_conditional_log_score = mean(data$log_score),
      failed_fits = sum(
        diagnostics$model_set == data$model_set[1] &
          diagnostics$prior == data$prior[1] &
          is.na(diagnostics$convergence)
      ),
      stringsAsFactors = FALSE
    )
  }))

  metrics <- metrics[
    order(
      match(metrics$model_set, expected_sets),
      match(metrics$prior, c("baseline", priors))
    ),
    , drop = FALSE
  ]
  rownames(metrics) <- NULL

  # Ranks compare predictor sets within the same prior family.
  rank_within_prior <- function(values, prior_values) {
    result <- rep(NA_real_, length(values))
    for (prior_name in priors) {
      rows <- prior_values == prior_name
      result[rows] <- rank(values[rows], ties.method = "min", na.last = "keep")
    }
    result
  }

  metrics$RMSE_rank_within_prior <- rank_within_prior(metrics$RMSE, metrics$prior)
  metrics$MAE_rank_within_prior <- rank_within_prior(metrics$MAE, metrics$prior)
  metrics$log_score_rank_within_prior <- rank_within_prior(
    -metrics$mean_conditional_log_score, metrics$prior
  )

  # Refit results should reproduce the existing full-five and baseline runs.
  compare_numeric <- function(x, y) {
    isTRUE(all.equal(
      unname(x), unname(y), tolerance = 1e-7,
      check.attributes = FALSE
    ))
  }

  original_predictions <- baseline_predictions
  original_predictions$Date <- as.Date(original_predictions$Date)
  original_diagnostics <- baseline_diagnostics
  original_diagnostics$Date <- as.Date(original_diagnostics$Date)

  reproduction_checks <- list()
  for (set_name in c("constant_gpd", "full_5")) {
    prior_names <- if (set_name == "constant_gpd") "baseline" else priors

    for (prior_name in prior_names) {
      original_model <- if (set_name == "constant_gpd") {
        "constant_gpd"
      } else {
        prior_name
      }

      new_rows <- predictions[
        predictions$model_set == set_name & predictions$prior == prior_name,
        , drop = FALSE
      ]
      old_rows <- original_predictions[
        original_predictions$model == original_model,
        , drop = FALSE
      ]
      new_rows <- new_rows[order(new_rows$Date), , drop = FALSE]
      old_rows <- old_rows[order(old_rows$Date), , drop = FALSE]

      same <- nrow(new_rows) == nrow(old_rows) &&
        nrow(new_rows) == length(forecast_dates)
      if (same) {
        same <- identical(new_rows$Date, old_rows$Date) &&
          compare_numeric(new_rows$observed_excess, old_rows$observed_excess) &&
          compare_numeric(new_rows$predicted_scale, old_rows$predicted_scale) &&
          compare_numeric(new_rows$predicted_mean_excess, old_rows$predicted_mean_excess) &&
          compare_numeric(new_rows$shape, old_rows$shape) &&
          compare_numeric(new_rows$log_score, old_rows$log_score)
      }

      reproduction_checks[[length(reproduction_checks) + 1L]] <- data.frame(
        check = paste("Matches original rolling pipeline:", set_name, prior_name),
        passed = isTRUE(same),
        detail = sprintf(
          "%d comparison forecasts; %d original forecasts",
          nrow(new_rows), nrow(old_rows)
        ),
        stringsAsFactors = FALSE
      )
    }
  }

  combination_grid <- expand.grid(
    model_set = expected_sets,
    prior = c("baseline", priors),
    stringsAsFactors = FALSE
  )
  combination_grid <- combination_grid[
    (combination_grid$model_set == "constant_gpd" &
       combination_grid$prior == "baseline") |
      (combination_grid$model_set != "constant_gpd" &
         combination_grid$prior != "baseline"),
    , drop = FALSE
  ]

  combination_counts <- vapply(seq_len(nrow(combination_grid)), function(i) {
    sum(
      predictions$model_set == combination_grid$model_set[i] &
        predictions$prior == combination_grid$prior[i]
    )
  }, integer(1))

  baseline_training <- original_diagnostics[
    original_diagnostics$model == "constant_gpd",
    c("Date", "n_train"), drop = FALSE
  ]
  baseline_training <- baseline_training[order(baseline_training$Date), ]
  comparison_training <- windows[
    order(windows$forecast_date),
    c("forecast_date", "n_training_exceedances"),
    drop = FALSE
  ]

  same_training_counts <-
    nrow(baseline_training) == length(forecast_dates) &&
    nrow(comparison_training) == length(forecast_dates) &&
    identical(baseline_training$Date, comparison_training$forecast_date) &&
    identical(
      as.integer(baseline_training$n_train),
      as.integer(comparison_training$n_training_exceedances)
    )

  all_fits_converged <-
    all(!is.na(diagnostics$convergence) & diagnostics$convergence == 0L) &&
    all(is.finite(predictions$predicted_scale) & predictions$predicted_scale > 0) &&
    all(is.finite(predictions$log_score))

  checks <- do.call(rbind, c(
    reproduction_checks,
    list(
      data.frame(
        check = "Every model was scored on every common forecast date",
        passed = nrow(predictions) == expected_rows &&
          all(combination_counts == length(forecast_dates)),
        detail = sprintf(
          "%d expected and observed rows; %d dates per model/prior combination",
          expected_rows, min(combination_counts)
        ),
        stringsAsFactors = FALSE
      ),
      data.frame(
        check = "Training counts match the original rolling folds",
        passed = same_training_counts,
        detail = sprintf("%d rolling folds compared", nrow(windows)),
        stringsAsFactors = FALSE
      ),
      data.frame(
        check = "All fits converged and predictions are finite and in support",
        passed = all_fits_converged,
        detail = sprintf(
          "%d converged fits; %d failed; %d non-finite log scores",
          sum(!is.na(diagnostics$convergence) & diagnostics$convergence == 0L),
          sum(is.na(diagnostics$convergence) | diagnostics$convergence != 0L),
          sum(!is.finite(predictions$log_score))
        ),
        stringsAsFactors = FALSE
      ),
      data.frame(
        check = "Training windows are strictly prior and meet minimum sample size",
        passed = all(windows$window_last_date < windows$forecast_date) &&
          all(windows$n_window_dates == config$rolling_observations) &&
          all(windows$n_training_exceedances >= config$minimum_training_exceedances),
        detail = sprintf(
          "%d folds; %d to %d training exceedances",
          nrow(windows), min(windows$n_training_exceedances),
          max(windows$n_training_exceedances)
        ),
        stringsAsFactors = FALSE
      )
    )
  ))

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  write_comparison_csv <- function(data, filename) {
    data_to_write <- data
    date_columns <- c(
      "Date", "forecast_date", "window_first_date", "window_last_date"
    )
    for (column in intersect(date_columns, names(data_to_write))) {
      data_to_write[[column]] <- format(
        as.Date(data_to_write[[column]]), "%Y-%m-%d"
      )
    }
    utils::write.csv(
      data_to_write,
      file.path(output_dir, filename),
      row.names = FALSE,
      na = ""
    )
  }

  write_comparison_csv(predictions, "reduced_model_comparison_predictions.csv")
  write_comparison_csv(metrics, "reduced_model_comparison_metrics.csv")
  write_comparison_csv(diagnostics, "reduced_model_comparison_fit_diagnostics.csv")
  write_comparison_csv(windows, "reduced_model_comparison_windows.csv")
  write_comparison_csv(checks, "reduced_model_comparison_validation.csv")

  saveRDS(
    list(
      model_sets = model_sets,
      predictions = predictions,
      metrics = metrics,
      diagnostics = diagnostics,
      windows = windows,
      validation = checks
    ),
    file.path(output_dir, "reduced_model_comparison.rds")
  )

  if (any(!checks$passed)) {
    failed <- checks$check[!checks$passed]
    stop(
      "Reduced-model comparison validation failed: ",
      paste(failed, collapse = "; ")
    )
  }

  list(
    model_sets = model_sets,
    predictions = predictions,
    metrics = metrics,
    diagnostics = diagnostics,
    windows = windows,
    validation = checks
  )
}
