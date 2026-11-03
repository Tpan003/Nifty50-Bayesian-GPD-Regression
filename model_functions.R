# Generalised Pareto likelihood, prior penalties, MAP estimation, and
# expanding/rolling out-of-sample evaluation for exceedance severity.

# Log density of a GPD observation z, with log(scale) supplied for numerical
# stability. The exponential distribution is the xi = 0 limiting case.
gpd_log_density <- function(z, log_sigma, xi) {
  if (length(z) != length(log_sigma) ||
      any(!is.finite(z)) || any(z <= 0) ||
      any(!is.finite(log_sigma)) || any(abs(log_sigma) > 25) ||
      length(xi) != 1L || !is.finite(xi)) {
    return(rep(-Inf, length(z)))
  }

  sigma <- exp(log_sigma)
  if (abs(xi) < 1e-7) {
    return(-log_sigma - z / sigma)
  }

  support <- 1 + xi * z / sigma
  if (any(!is.finite(support)) || any(support <= 0)) {
    return(rep(-Inf, length(z)))
  }

  -log_sigma - (1 / xi + 1) * log1p(xi * z / sigma)
}

# All prior hyperparameters are defined here and saved with the model outputs.
gpd_prior_specification <- function() {
  list(
    beta_cauchy_scale = 1,
    beta_lasso_rate = 1,
    beta_lasso_epsilon = 1e-3,
    beta_ridge_precision = 1,
    g_prior_g = "n",
    xi_cauchy_scale = 1,
    xi_lower = -5,
    xi_upper = 1 - 1e-6
  )
}

shape_prior_penalty <- function(xi, scale = 1) {
  log1p((xi / scale)^2)
}

# Slope priors are applied to standardised predictors. The intercept is not
# penalised. The lasso penalty is smoothed slightly for stable optimisation.
prior_penalty <- function(beta, prior, x, prior_spec = gpd_prior_specification()) {
  slopes <- beta[-1L]

  if (prior == "cauchy") {
    scale <- prior_spec$beta_cauchy_scale
    return(sum(log1p((slopes / scale)^2)))
  }

  if (prior == "lasso") {
    epsilon <- prior_spec$beta_lasso_epsilon
    return(
      prior_spec$beta_lasso_rate *
        sum(sqrt(slopes^2 + epsilon^2) - epsilon)
    )
  }

  if (prior == "ridge") {
    return(prior_spec$beta_ridge_precision * sum(slopes^2) / 2)
  }

  if (prior == "g-prior") {
    if (length(slopes) == 0L) return(0)
    design <- x[, -1L, drop = FALSE]
    g <- if (identical(prior_spec$g_prior_g, "n")) {
      nrow(design)
    } else {
      as.numeric(prior_spec$g_prior_g)
    }
    if (!is.finite(g) || g <= 0) stop("The g-prior scale must be positive.")
    return(as.numeric(crossprod(slopes, crossprod(design) %*% slopes)) / (2 * g))
  }

  stop("Unknown prior: ", prior)
}

gpd_map_objective <- function(p, z, x, prior = NULL,
                              prior_spec = gpd_prior_specification()) {
  beta <- p[seq_len(ncol(x))]
  xi <- p[length(p)]
  log_sigma <- drop(x %*% beta)
  log_likelihood <- gpd_log_density(z, log_sigma, xi)

  if (any(!is.finite(log_likelihood))) return(1e100)

  penalty <- if (is.null(prior)) {
    0
  } else {
    prior_penalty(beta, prior, x, prior_spec)
  }

  -sum(log_likelihood) + penalty +
    shape_prior_penalty(xi, prior_spec$xi_cauchy_scale)
}

# Analytic gradient of the negative log-posterior objective.
gpd_map_gradient <- function(p, z, x, prior = NULL,
                             prior_spec = gpd_prior_specification()) {
  beta <- p[seq_len(ncol(x))]
  xi <- p[length(p)]
  eta <- drop(x %*% beta)
  u <- z * exp(-eta)
  support <- 1 + xi * u

  if (any(!is.finite(support)) || any(support <= 0)) {
    return(rep(0, length(p)))
  }

  d_eta <- -1 + (1 + xi) * u / support
  d_xi <- if (abs(xi) < 1e-4) {
    -u + 0.5 * u^2 + xi * (u^2 - (2 / 3) * u^3)
  } else {
    log1p(xi * u) / xi^2 - (1 + 1 / xi) * u / support
  }

  grad_beta <- -drop(crossprod(x, d_eta))
  slopes <- beta[-1L]

  if (!is.null(prior) && prior == "cauchy") {
    scale <- prior_spec$beta_cauchy_scale
    grad_beta[-1L] <- grad_beta[-1L] +
      2 * slopes / (scale^2 + slopes^2)
  } else if (!is.null(prior) && prior == "lasso") {
    epsilon <- prior_spec$beta_lasso_epsilon
    grad_beta[-1L] <- grad_beta[-1L] +
      prior_spec$beta_lasso_rate * slopes / sqrt(slopes^2 + epsilon^2)
  } else if (!is.null(prior) && prior == "ridge") {
    grad_beta[-1L] <- grad_beta[-1L] +
      prior_spec$beta_ridge_precision * slopes
  } else if (!is.null(prior) && prior == "g-prior" && length(slopes) > 0L) {
    design <- x[, -1L, drop = FALSE]
    g <- if (identical(prior_spec$g_prior_g, "n")) {
      nrow(design)
    } else {
      as.numeric(prior_spec$g_prior_g)
    }
    if (!is.finite(g) || g <= 0) stop("The g-prior scale must be positive.")
    grad_beta[-1L] <- grad_beta[-1L] +
      drop(crossprod(design) %*% slopes) / g
  }

  xi_scale <- prior_spec$xi_cauchy_scale
  c(
    grad_beta,
    -sum(d_xi) + 2 * xi / (xi_scale^2 + xi^2)
  )
}

fit_constant_gpd <- function(z, prior_spec = gpd_prior_specification()) {
  if (length(z) == 0L || any(!is.finite(z)) || any(z <= 0)) {
    stop("Training excesses must be positive and finite.")
  }

  x <- matrix(1, nrow = length(z), ncol = 1L)
  objective <- function(p) {
    gpd_map_objective(p, z, x, prior = NULL, prior_spec = prior_spec)
  }
  gradient <- function(p) {
    gpd_map_gradient(p, z, x, prior = NULL, prior_spec = prior_spec)
  }

  start <- c(log(mean(z)), 0.05)
  fit <- stats::nlminb(
    start = start,
    objective = objective,
    gradient = gradient,
    lower = c(-10, prior_spec$xi_lower),
    upper = c(10, prior_spec$xi_upper),
    control = list(iter.max = 2000, eval.max = 4000,
                   rel.tol = 1e-10, x.tol = 1e-8)
  )

  start_value <- objective(start)
  if (fit$convergence != 0L || !is.finite(fit$objective) ||
      fit$objective >= 1e99 ||
      fit$objective > start_value + 1e-7 * max(1, abs(start_value))) {
    stop("Constant GPD optimisation failed: ", fit$message)
  }

  list(
    log_sigma = fit$par[1L],
    xi = fit$par[2L],
    convergence = fit$convergence,
    objective = fit$objective,
    iterations = fit$iterations,
    evaluations = fit$evaluations,
    message = fit$message
  )
}

fit_gpd_map <- function(z, predictors, prior = "ridge",
                        prior_spec = gpd_prior_specification()) {
  predictors <- as.matrix(predictors)
  storage.mode(predictors) <- "double"

  if (nrow(predictors) != length(z) || nrow(predictors) < 2L ||
      ncol(predictors) < 1L || any(!is.finite(predictors)) ||
      any(!is.finite(z)) || any(z <= 0)) {
    stop("Invalid GPD training data.")
  }

  allowed_priors <- c("cauchy", "lasso", "ridge", "g-prior")
  if (length(prior) != 1L || !prior %in% allowed_priors) {
    stop("Unknown prior: ", prior)
  }

  predictor_means <- colMeans(predictors)
  predictor_sds <- apply(predictors, 2L, stats::sd)
  if (any(!is.finite(predictor_sds)) || any(predictor_sds <= 1e-10)) {
    stop("At least one training predictor has zero or near-zero variation.")
  }

  x <- cbind(
    Intercept = 1,
    scale(predictors, center = predictor_means, scale = predictor_sds)
  )
  constant_fit <- fit_constant_gpd(z, prior_spec)
  start <- c(constant_fit$log_sigma, rep(0, ncol(x) - 1L), constant_fit$xi)

  objective <- function(p) {
    gpd_map_objective(p, z, x, prior, prior_spec)
  }
  gradient <- function(p) {
    gpd_map_gradient(p, z, x, prior, prior_spec)
  }

  fit <- stats::nlminb(
    start = start,
    objective = objective,
    gradient = gradient,
    lower = c(rep(-15, ncol(x)), prior_spec$xi_lower),
    upper = c(rep(15, ncol(x)), prior_spec$xi_upper),
    control = list(iter.max = 2000, eval.max = 4000,
                   rel.tol = 1e-10, x.tol = 1e-8)
  )

  start_value <- objective(start)
  if (fit$convergence != 0L || !is.finite(fit$objective) ||
      fit$objective >= 1e99 ||
      fit$objective > start_value + 1e-7 * max(1, abs(start_value))) {
    stop("MAP optimisation failed: ", fit$message)
  }

  beta <- fit$par[seq_len(ncol(x))]
  xi <- fit$par[length(fit$par)]
  fitted_log_scale <- drop(x %*% beta)
  if (any(!is.finite(gpd_log_density(z, fitted_log_scale, xi)))) {
    stop("The fitted parameters are outside the GPD support.")
  }

  structure(
    list(
      beta = beta,
      xi = xi,
      prior = prior,
      means = predictor_means,
      sds = predictor_sds,
      predictor_names = colnames(predictors),
      convergence = fit$convergence,
      objective = fit$objective,
      iterations = fit$iterations,
      evaluations = fit$evaluations,
      message = fit$message,
      prior_specification = prior_spec,
      n = length(z)
    ),
    class = "gpd_map"
  )
}

predict.gpd_map <- function(object, newdata, ...) {
  newdata <- as.matrix(newdata)
  storage.mode(newdata) <- "double"

  if (ncol(newdata) != length(object$means) || any(!is.finite(newdata))) {
    stop("Prediction data has invalid values or the wrong number of columns.")
  }
  if (!is.null(object$predictor_names) && !is.null(colnames(newdata)) &&
      !identical(colnames(newdata), object$predictor_names)) {
    stop("Prediction columns must match the training columns and their order.")
  }

  x <- cbind(
    Intercept = 1,
    scale(newdata, center = object$means, scale = object$sds)
  )
  exp(drop(x %*% object$beta))
}

rolling_evaluation <- function(merged_data, config) {
  predictors <- c(
    "nifty_vol_lag1", "sp500_vol_lag1", "gold_vol_lag1",
    "india_vix_lag1", "usdinr_vol_lag1"
  )
  required_columns <- c("Date", "exceedance", "excess", predictors)
  if (!all(required_columns %in% names(merged_data))) {
    stop("Merged data is missing columns required for rolling evaluation.")
  }
  if (config$rolling_observations < 1L ||
      config$minimum_training_exceedances < 1L ||
      config$forecast_step < 1L) {
    stop("Rolling configuration values must be positive integers.")
  }

  prior_spec <- config$prior_specification
  dates <- sort(unique(as.Date(merged_data$Date)))
  first_origin <- min(config$rolling_observations + 1L, length(dates))
  forecast_dates <- dates[seq.int(
    first_origin, length(dates), by = config$forecast_step
  )]

  predictions <- list()
  diagnostics <- list()
  model_names <- c("constant_gpd", "cauchy", "lasso", "ridge", "g-prior")

  add_diagnostic <- function(date, model, n_train, fit = NULL, error = NULL,
                             max_slope = NA_real_) {
    failed <- !is.null(error)
    evaluations <- if (!failed && !is.null(fit$evaluations)) {
      fit$evaluations
    } else {
      c("function" = NA_integer_, gradient = NA_integer_)
    }

    diagnostics[[length(diagnostics) + 1L]] <<- data.frame(
      Date = as.Date(date),
      model = model,
      n_train = as.integer(n_train),
      convergence = if (failed) NA_integer_ else as.integer(fit$convergence),
      objective = if (failed) NA_real_ else as.numeric(fit$objective),
      iterations = if (failed) NA_integer_ else as.integer(fit$iterations),
      function_evaluations = unname(evaluations["function"]),
      gradient_evaluations = unname(evaluations["gradient"]),
      max_abs_slope = as.numeric(max_slope),
      message = if (failed) conditionMessage(error) else as.character(fit$message),
      stringsAsFactors = FALSE
    )
  }

  for (date in forecast_dates) {
    previous_dates <- sort(unique(merged_data$Date[merged_data$Date < date]))
    window_dates <- tail(previous_dates, config$rolling_observations)
    training_window <- merged_data[
      merged_data$Date %in% window_dates, , drop = FALSE
    ]

    training_events <- training_window[
      training_window$exceedance == 1L &
        is.finite(training_window$excess) &
        training_window$excess > 0 &
        complete.cases(training_window[, predictors, drop = FALSE]),
      , drop = FALSE
    ]

    if (nrow(training_events) < config$minimum_training_exceedances) next

    target <- merged_data[
      merged_data$Date == date &
        merged_data$exceedance == 1L &
        is.finite(merged_data$excess) &
        merged_data$excess > 0 &
        complete.cases(merged_data[, predictors, drop = FALSE]),
      , drop = FALSE
    ]
    if (nrow(target) == 0L) next

    baseline <- tryCatch(
      fit_constant_gpd(training_events$excess, prior_spec),
      error = function(e) e
    )

    if (inherits(baseline, "error")) {
      add_diagnostic(date, "constant_gpd", nrow(training_events), error = baseline,
                     max_slope = 0)
      for (model in model_names[-1L]) {
        add_diagnostic(
          date, model, nrow(training_events),
          error = simpleError("Skipped because constant GPD initialisation failed.")
        )
      }
      next
    }

    add_diagnostic(date, "constant_gpd", nrow(training_events),
                   fit = baseline, max_slope = 0)

    baseline_scale <- exp(baseline$log_sigma)
    baseline_mean <- if (baseline$xi < 1) {
      baseline_scale / (1 - baseline$xi)
    } else {
      NA_real_
    }

    predictions[[length(predictions) + 1L]] <- data.frame(
      Date = target$Date,
      observed_excess = target$excess,
      model = "constant_gpd",
      predicted_scale = baseline_scale,
      predicted_mean_excess = baseline_mean,
      shape = baseline$xi,
      log_score = gpd_log_density(
        target$excess, rep(baseline$log_sigma, nrow(target)), baseline$xi
      ),
      stringsAsFactors = FALSE
    )

    for (prior in model_names[-1L]) {
      fitted <- tryCatch(
        fit_gpd_map(
          training_events$excess,
          training_events[, predictors, drop = FALSE],
          prior = prior,
          prior_spec = prior_spec
        ),
        error = function(e) e
      )

      if (inherits(fitted, "error")) {
        add_diagnostic(date, prior, nrow(training_events), error = fitted)
        next
      }

      predicted_scale <- predict(
        fitted, target[, predictors, drop = FALSE]
      )
      predicted_mean <- if (fitted$xi < 1) {
        predicted_scale / (1 - fitted$xi)
      } else {
        rep(NA_real_, length(predicted_scale))
      }

      predictions[[length(predictions) + 1L]] <- data.frame(
        Date = target$Date,
        observed_excess = target$excess,
        model = prior,
        predicted_scale = predicted_scale,
        predicted_mean_excess = predicted_mean,
        shape = fitted$xi,
        log_score = gpd_log_density(
          target$excess, log(predicted_scale), fitted$xi
        ),
        stringsAsFactors = FALSE
      )

      add_diagnostic(
        date, prior, nrow(training_events), fit = fitted,
        max_slope = max(abs(fitted$beta[-1L]))
      )
    }
  }

  empty_predictions <- data.frame(
    Date = as.Date(character()),
    observed_excess = numeric(),
    model = character(),
    predicted_scale = numeric(),
    predicted_mean_excess = numeric(),
    shape = numeric(),
    log_score = numeric()
  )

  result <- if (length(predictions) > 0L) {
    do.call(rbind, predictions)
  } else {
    empty_predictions
  }
  diagnostic_table <- if (length(diagnostics) > 0L) {
    do.call(rbind, diagnostics)
  } else {
    data.frame()
  }

  attr(result, "fit_diagnostics") <- diagnostic_table
  result
}

score_predictions <- function(predictions) {
  if (nrow(predictions) == 0L) return(data.frame())

  by_model <- split(predictions, predictions$model)
  do.call(rbind, lapply(names(by_model), function(model_name) {
    data <- by_model[[model_name]]
    valid_mean <- is.finite(data$predicted_mean_excess)

    data.frame(
      model = model_name,
      forecasts = nrow(data),
      forecasts_with_finite_mean = sum(valid_mean),
      nonfinite_log_scores = sum(!is.finite(data$log_score)),
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
      stringsAsFactors = FALSE
    )
  }))
}
