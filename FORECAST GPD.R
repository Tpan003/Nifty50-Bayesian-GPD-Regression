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

# ==============================================================================
# STEP 0: PACKAGE VERIFICATION
# ==============================================================================
required <- c("readxl", "ggplot2", "knitr", "rmarkdown")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))

# ==============================================================================
# STEP 1: DATA IMPORT, VALIDATION, DATE ALIGNMENT & EXCEEDANCE CONSTRUCTION
# Source file: R/data_functions.R
# ==============================================================================
# Data import, validation, date alignment and exceedance construction.

read_market_sheets <- function(path) {
  sheets <- readxl::excel_sheets(path)
  expected <- c("NIFTY 50 DATA INR", "S&P500 DATA USD", "GOLD DATA USD",
                "INDIA VIX DATA INR", "USD-INR DATA")
  absent <- setdiff(expected, sheets)
  if (length(absent)) stop("Workbook is missing sheets: ", paste(absent, collapse = ", "))
  # The workbook uses a period in some numeric cells to mean missing. Pass
  # this sentinel to readxl and suppress its expected per-cell coercion notices;
  # downstream date, response and predictor checks still validate the imports.
  read_sheet <- function(s) suppressMessages(suppressWarnings(
    readxl::read_excel(path, sheet = s, na = c("", "."))))
  stats::setNames(lapply(expected, read_sheet),
                  c("nifty", "sp500", "gold", "vix", "usdinr"))
}

check_dates <- function(x, label) {
  if (!"Date" %in% names(x)) stop(label, " has no Date column")
  x$Date <- as.Date(x$Date)
  if (anyNA(x$Date)) stop(label, " contains missing or unparseable dates")
  if (anyDuplicated(x$Date)) stop(label, " contains duplicate dates")
  x[order(x$Date), , drop = FALSE]
}

# At forecast date t, use the most recent source observation strictly before t.
# This handles holidays in each market without requiring identical calendars.
latest_prior <- function(target_dates, source_dates, values) {
  stopifnot(length(source_dates) == length(values))
  ok <- !is.na(source_dates) & is.finite(values)
  source_dates <- as.Date(source_dates[ok]); values <- values[ok]
  o <- order(source_dates); source_dates <- source_dates[o]; values <- values[o]
  idx <- findInterval(as.numeric(as.Date(target_dates) - 1), as.numeric(source_dates))
  ans <- rep(NA_real_, length(target_dates)); have <- idx > 0L
  ans[have] <- values[idx[have]]
  ans
}

prepare_and_merge <- function(path, threshold = 2) {
  stopifnot(length(threshold) == 1L, is.finite(threshold), threshold > 0)
  s <- read_market_sheets(path)
  n <- s$nifty
  required <- c("DATE", "Log return (%)", "Daily empirical volatility (%)")
  if (!all(required %in% names(n))) stop("Nifty sheet is missing required columns")
  n <- data.frame(Date = as.Date(n$DATE),
                  nifty_return = as.numeric(n[["Log return (%)"]]),
                  nifty_vol = as.numeric(n[["Daily empirical volatility (%)"]]))
  n <- check_dates(n, "Nifty sheet")
  n$nifty_vol_lag1 <- latest_prior(n$Date, n$Date, n$nifty_vol)
  sources <- list(
    sp500 = list(sheet = s$sp500, column = "Daily empirical volatility (%)", name = "sp500_vol"),
    gold = list(sheet = s$gold, column = "Daily empirical volatility (%)", name = "gold_vol"),
    vix = list(sheet = s$vix, column = "India VIX level (points)", name = "india_vix"),
    usdinr = list(sheet = s$usdinr, column = "Daily empirical volatility (%)", name = "usdinr_vol")
  )
  for (key in names(sources)) {
    z <- sources[[key]]; dcol <- if ("DATE" %in% names(z$sheet)) "DATE" else "Date"
    if (!z$column %in% names(z$sheet) || !dcol %in% names(z$sheet)) stop("Missing columns on ", key, " sheet")
    source <- check_dates(data.frame(Date = as.Date(z$sheet[[dcol]]),
                                     value = as.numeric(z$sheet[[z$column]])), key)
    n[[paste0(z$name, "_lag1")]] <- latest_prior(n$Date, source$Date, source$value)
  }
  n$threshold <- threshold
  n$exceedance <- as.integer(is.finite(n$nifty_return) & n$nifty_return < -threshold)
  n$excess <- ifelse(n$exceedance == 1L, -n$nifty_return - threshold, NA_real_)
  # Cross-check stored workbook excesses when present; computed response is authoritative.
  original_name <- paste0("Excess loss beyond ", threshold, "% threshold (pp)")
  if (original_name %in% names(s$nifty)) {
    stored <- data.frame(Date = as.Date(s$nifty$DATE),
                         stored_excess = as.numeric(s$nifty[[original_name]]))
    stored <- check_dates(stored, "stored Nifty excess column")
    stored_by_date <- stored$stored_excess[match(n$Date, stored$Date)]
    comparable <- is.finite(stored_by_date) & is.finite(n$excess)
    if (any(comparable) && max(abs(stored_by_date[comparable] - n$excess[comparable])) > 1e-6)
      warning("Stored excess-loss column differs from recomputed values; recomputed values are used.")
  }
  n
}

build_extremes <- function(merged) {
  predictors <- c("nifty_vol", "sp500_vol", "gold_vol", "india_vix", "usdinr_vol")
  names(predictors) <- paste0(predictors, "_lag1")
  x <- merged[merged$exceedance == 1L & is.finite(merged$excess) & merged$excess > 0, , drop = FALSE]
  for (nm in names(predictors)) x[[nm]] <- x[[predictors[[nm]]]]
  x <- x[order(x$Date), , drop = FALSE]
  x
}

# ==============================================================================
# STEP 2: GPD LIKELIHOOD, PRIOR PENALTIES, MAP FITTING & ROLLING EVALUATION
# Source file: R/model_functions.R
# ==============================================================================
# GPD likelihood, prior-regularized MAP estimators and rolling evaluation.

gpd_log_density <- function(z, log_sigma, xi) {
  if (length(z) != length(log_sigma) || any(!is.finite(z)) || any(z <= 0) ||
      any(!is.finite(log_sigma)) || any(abs(log_sigma) > 25) ||
      length(xi) != 1L || !is.finite(xi))
    return(rep(-Inf, length(z)))
  sigma <- exp(log_sigma)
  if (abs(xi) < 1e-7) return(-log_sigma - z / sigma)
  support <- 1 + xi * z / sigma
  if (any(!is.finite(support)) || any(support <= 0)) return(rep(-Inf, length(z)))
  -log_sigma - (1 / xi + 1) * log1p(xi * z / sigma)
}

# Prior hyperparameters are returned centrally and saved with the outputs.
gpd_prior_specification <- function() {
  list(beta_cauchy_scale = 1, beta_lasso_rate = 1,
       beta_lasso_epsilon = 1e-3, beta_ridge_precision = 1,
       g_prior_g = "n", xi_cauchy_scale = 1,
       xi_lower = -5, xi_upper = 1 - 1e-6)
}

# Negative log density for the paper's truncated Cauchy prior on xi.
# The normalizing constant is independent of xi and can be omitted for MAP.
shape_prior_penalty <- function(xi, scale = 1) log1p((xi / scale)^2)

# Prior penalties apply to standardized slopes; the intercept is unpenalized.
prior_penalty <- function(beta, prior, X, prior_spec = gpd_prior_specification()) {
  b <- beta[-1L]
  if (prior == "cauchy") return(sum(log1p((b / prior_spec$beta_cauchy_scale)^2)))
  # [FIX B-09] use the same smoothed L1 penalty that the optimizer actually uses
  if (prior == "lasso") return(prior_spec$beta_lasso_rate *
                                 sum(sqrt(b^2 + prior_spec$beta_lasso_epsilon^2) - prior_spec$beta_lasso_epsilon))
  if (prior == "ridge") return(prior_spec$beta_ridge_precision * sum(b^2) / 2)
  if (prior == "g-prior") {
    if (!length(b)) return(0)
    design <- X[, -1L, drop = FALSE]
    gram <- crossprod(design)
    g <- if (identical(prior_spec$g_prior_g, "n")) nrow(design) else as.numeric(prior_spec$g_prior_g)
    if (!is.finite(g) || g <= 0) stop("g-prior scale must be positive")
    return(as.numeric(crossprod(b, gram %*% b)) / (2 * g))
  }
  stop("Unknown prior")
}

gpd_map_objective <- function(p, z, X, prior = NULL,
                              prior_spec = gpd_prior_specification()) {
  beta <- p[seq_len(ncol(X))]
  xi <- p[length(p)]
  ll <- gpd_log_density(z, drop(X %*% beta), xi)
  if (any(!is.finite(ll))) return(1e100)
  penalty <- if (is.null(prior)) 0 else prior_penalty(beta, prior, X, prior_spec)
  -sum(ll) + penalty + shape_prior_penalty(xi, prior_spec$xi_cauchy_scale)
}

# Analytic gradient for the smooth GPD likelihood and priors. The lasso's
# absolute-value penalty is smoothed at epsilon=1e-3 for stable optimization;
# this differs from the exact L1 penalty by at most p * epsilon.
gpd_map_gradient <- function(p, z, X, prior = NULL,
                             prior_spec = gpd_prior_specification()) {
  beta <- p[seq_len(ncol(X))]
  xi <- p[length(p)]
  eta <- drop(X %*% beta)
  u <- z * exp(-eta)
  support <- 1 + xi * u
  if (any(!is.finite(support)) || any(support <= 0)) return(rep(0, length(p)))
  d_eta <- -1 + (1 + xi) * u / support
  d_xi <- if (abs(xi) < 1e-4) {
    -u + 0.5 * u^2 + xi * (u^2 - (2 / 3) * u^3)
  } else {
    log1p(xi * u) / xi^2 - (1 + 1 / xi) * u / support
  }
  grad_beta <- -drop(crossprod(X, d_eta))
  slopes <- beta[-1L]
  if (!is.null(prior) && prior == "cauchy") {
    scale <- prior_spec$beta_cauchy_scale
    grad_beta[-1L] <- grad_beta[-1L] + 2 * slopes / (scale^2 + slopes^2)
  } else if (!is.null(prior) && prior == "lasso") {
    epsilon <- prior_spec$beta_lasso_epsilon
    grad_beta[-1L] <- grad_beta[-1L] + prior_spec$beta_lasso_rate * slopes / sqrt(slopes^2 + epsilon^2)
  } else if (!is.null(prior) && prior == "ridge") {
    grad_beta[-1L] <- grad_beta[-1L] + prior_spec$beta_ridge_precision * slopes
  } else if (!is.null(prior) && prior == "g-prior" && length(slopes)) {
    design <- X[, -1L, drop = FALSE]
    g <- if (identical(prior_spec$g_prior_g, "n")) nrow(design) else as.numeric(prior_spec$g_prior_g)
    if (!is.finite(g) || g <= 0) stop("g-prior scale must be positive")
    grad_beta[-1L] <- grad_beta[-1L] + drop(crossprod(design) %*% slopes) / g
  }
  xi_scale <- prior_spec$xi_cauchy_scale
  c(grad_beta, -sum(d_xi) + 2 * xi / (xi_scale^2 + xi^2))
}

fit_constant_gpd <- function(z, prior_spec = gpd_prior_specification()) {
  if (!length(z) || any(!is.finite(z)) || any(z <= 0)) stop("Invalid GPD training excesses")
  design <- matrix(1, nrow = length(z), ncol = 1L)
  objective <- function(p) gpd_map_objective(p, z, design, prior_spec = prior_spec)
  gradient <- function(p) gpd_map_gradient(p, z, design, prior_spec = prior_spec)
  start <- c(log(mean(z)), 0.05)
  fit <- stats::nlminb(start, objective, gradient,
                       lower = c(-10, prior_spec$xi_lower),
                       upper = c(10, prior_spec$xi_upper),
                       control = list(iter.max = 2000, eval.max = 4000,
                                      rel.tol = 1e-10, x.tol = 1e-8))
  if (fit$convergence != 0L || !is.finite(fit$objective) || fit$objective >= 1e99 ||
      fit$objective > objective(start) + 1e-7 * max(1, abs(objective(start))))
    stop("Constant GPD optimizer failed to improve a feasible starting fit: ", fit$message)
  list(log_sigma = fit$par[1L], xi = fit$par[2L], convergence = fit$convergence,
       objective = fit$objective, iterations = fit$iterations,
       evaluations = fit$evaluations, message = fit$message)
}

fit_gpd_map <- function(z, predictors, prior = "ridge",
                        prior_spec = gpd_prior_specification()) {
  predictors <- as.matrix(predictors)
  storage.mode(predictors) <- "double"
  if (nrow(predictors) != length(z) || nrow(predictors) < 2L || ncol(predictors) < 1L ||
      any(!is.finite(predictors)) || any(!is.finite(z)) || any(z <= 0))
    stop("Invalid GPD training data")
  priors <- c("cauchy", "lasso", "ridge", "g-prior")
  if (length(prior) != 1L || !prior %in% priors) stop("Unknown prior")
  means <- colMeans(predictors)
  sds <- apply(predictors, 2L, stats::sd)
  if (any(!is.finite(sds)) || any(sds <= 1e-10)) stop("Training predictor has zero variation")
  X <- cbind(Intercept = rep(1, nrow(predictors)),
             scale(predictors, center = means, scale = sds))
  constant <- fit_constant_gpd(z, prior_spec)
  start <- c(constant$log_sigma, rep(0, ncol(X) - 1L), constant$xi)
  objective <- function(p) gpd_map_objective(p, z, X, prior, prior_spec)
  gradient <- function(p) gpd_map_gradient(p, z, X, prior, prior_spec)
  fit <- stats::nlminb(start, objective, gradient,
                       lower = c(rep(-15, ncol(X)), prior_spec$xi_lower),
                       upper = c(rep(15, ncol(X)), prior_spec$xi_upper),
                       control = list(iter.max = 2000, eval.max = 4000,
                                      rel.tol = 1e-10, x.tol = 1e-8))
  initial_value <- objective(start)
  if (fit$convergence != 0L || !is.finite(fit$objective) || fit$objective >= 1e99 ||
      fit$objective > initial_value + 1e-7 * max(1, abs(initial_value)))
    stop("MAP optimizer failed to improve a feasible starting fit: ", fit$message)
  beta <- fit$par[seq_len(ncol(X))]
  xi <- fit$par[length(fit$par)]
  log_sigma <- drop(X %*% beta)
  if (any(!is.finite(gpd_log_density(z, log_sigma, xi))))
    stop("MAP optimizer returned parameters outside GPD support")
  structure(list(beta = beta, xi = xi, prior = prior, means = means, sds = sds,
                 predictor_names = colnames(predictors), convergence = fit$convergence,
                 objective = fit$objective, iterations = fit$iterations,
                 evaluations = fit$evaluations, message = fit$message,
                 prior_specification = prior_spec, n = length(z)), class = "gpd_map")
}

predict.gpd_map <- function(object, newdata, ...) {
  x <- as.matrix(newdata)
  storage.mode(x) <- "double"
  if (ncol(x) != length(object$means) || any(!is.finite(x)))
    stop("Invalid prediction predictors")
  if (!is.null(object$predictor_names) && !is.null(colnames(x)) &&
      !identical(colnames(x), object$predictor_names))
    stop("Prediction columns differ from training columns or order")
  X <- cbind(Intercept = rep(1, nrow(x)),
             scale(x, center = object$means, scale = object$sds))
  exp(drop(X %*% object$beta))
}

rolling_evaluation <- function(merged, cfg) {
  predictors <- c("nifty_vol_lag1", "sp500_vol_lag1", "gold_vol_lag1",
                  "india_vix_lag1", "usdinr_vol_lag1")
  required <- c("Date", "exceedance", "excess", predictors)
  if (!all(required %in% names(merged))) stop("Merged data is missing rolling-model columns")
  if (cfg$rolling_observations < 1L || cfg$minimum_training_exceedances < 1L ||
      cfg$forecast_step < 1L) stop("Rolling configuration values must be positive integers")
  prior_spec <- if (is.null(cfg$prior_specification)) gpd_prior_specification() else cfg$prior_specification
  dates <- sort(unique(merged$Date))
  first_origin <- min(cfg$rolling_observations + 1L, length(dates))
  dates <- dates[seq.int(first_origin, length(dates), by = cfg$forecast_step)]
  predictions <- list()
  diagnostics <- list()
  model_names <- c("constant_gpd", "cauchy", "lasso", "ridge", "g-prior")

  add_diagnostic <- function(date, model, n_train, fit = NULL, error = NULL, max_slope = NA_real_) {
    failed <- !is.null(error)
    evaluations <- if (!failed && !is.null(fit$evaluations)) fit$evaluations else c("function" = NA_integer_, gradient = NA_integer_)
    diagnostics[[length(diagnostics) + 1L]] <<- data.frame(
      Date = as.Date(date), model = model, n_train = as.integer(n_train),
      convergence = if (failed) NA_integer_ else as.integer(fit$convergence),
      objective = if (failed) NA_real_ else as.numeric(fit$objective),
      iterations = if (failed) NA_integer_ else as.integer(fit$iterations),
      function_evaluations = unname(evaluations["function"]),
      gradient_evaluations = unname(evaluations["gradient"]),
      max_abs_slope = as.numeric(max_slope),
      message = if (failed) conditionMessage(error) else as.character(fit$message),
      stringsAsFactors = FALSE)
  }

  # [FIX B-06] iterate by index so `date` keeps its Date class (as.Date(<numeric>)
  # without an origin fails on R < 4.3.0).
  for (di in seq_along(dates)) {
    date <- dates[di]
    earlier <- sort(unique(merged$Date[merged$Date < date]))
    keep_dates <- tail(earlier, cfg$rolling_observations)
    train <- merged[merged$Date %in% keep_dates, , drop = FALSE]
    train_ext <- train[train$exceedance == 1L & is.finite(train$excess) & train$excess > 0 &
                         complete.cases(train[, predictors, drop = FALSE]), , drop = FALSE]
    if (nrow(train_ext) < cfg$minimum_training_exceedances) next
    target <- merged[merged$Date == date & merged$exceedance == 1L &
                       is.finite(merged$excess) & merged$excess > 0 &
                       complete.cases(merged[, predictors, drop = FALSE]), , drop = FALSE]
    if (!nrow(target)) next

    base <- tryCatch(fit_constant_gpd(train_ext$excess, prior_spec), error = function(e) e)
    if (inherits(base, "error")) {
      add_diagnostic(date, "constant_gpd", nrow(train_ext), error = base, max_slope = 0)
      for (pr in model_names[-1L])
        add_diagnostic(date, pr, nrow(train_ext), error = simpleError("Skipped: constant GPD initialization failed"))
      next
    }
    add_diagnostic(date, "constant_gpd", nrow(train_ext), fit = base, max_slope = 0)
    base_ll <- gpd_log_density(target$excess, rep(base$log_sigma, nrow(target)), base$xi)
    base_scale <- exp(base$log_sigma)
    base_mean <- if (base$xi < 1) base_scale / (1 - base$xi) else rep(NA_real_, nrow(target))
    predictions[[length(predictions) + 1L]] <- data.frame(
      Date = target$Date, observed_excess = target$excess, model = "constant_gpd",
      predicted_scale = base_scale, predicted_mean_excess = base_mean,
      shape = base$xi, log_score = base_ll)

    for (pr in model_names[-1L]) {
      fitted <- tryCatch(fit_gpd_map(train_ext$excess, train_ext[, predictors, drop = FALSE],
                                     pr, prior_spec),
                         error = function(e) e)
      if (inherits(fitted, "error")) {
        add_diagnostic(date, pr, nrow(train_ext), error = fitted)
        next
      }
      scale_hat <- predict(fitted, target[, predictors, drop = FALSE])
      ll <- gpd_log_density(target$excess, log(scale_hat), fitted$xi)
      mean_hat <- if (fitted$xi < 1) scale_hat / (1 - fitted$xi) else rep(NA_real_, length(scale_hat))
      predictions[[length(predictions) + 1L]] <- data.frame(
        Date = target$Date, observed_excess = target$excess, model = pr,
        predicted_scale = scale_hat, predicted_mean_excess = mean_hat,
        shape = fitted$xi, log_score = ll)
      add_diagnostic(date, pr, nrow(train_ext), fit = fitted,
                     max_slope = max(abs(fitted$beta[-1L])))
    }
  }
  empty_predictions <- data.frame(Date = as.Date(character()), observed_excess = numeric(),
                                 model = character(), predicted_scale = numeric(),
                                 predicted_mean_excess = numeric(), shape = numeric(), log_score = numeric())
  result <- if (length(predictions)) do.call(rbind, predictions) else empty_predictions
  diag_table <- if (length(diagnostics)) do.call(rbind, diagnostics) else data.frame()
  attr(result, "fit_diagnostics") <- diag_table
  result
}

score_predictions <- function(predictions) {
  if (!nrow(predictions)) return(data.frame())
  split_data <- split(predictions, predictions$model)
  do.call(rbind, lapply(names(split_data), function(nm) {
    d <- split_data[[nm]]
    valid_mean <- is.finite(d$predicted_mean_excess)
    data.frame(model = nm, forecasts = nrow(d), forecasts_with_finite_mean = sum(valid_mean),
               nonfinite_log_scores = sum(!is.finite(d$log_score)),
               RMSE = if (any(valid_mean)) sqrt(mean((d$observed_excess[valid_mean] - d$predicted_mean_excess[valid_mean])^2)) else NA_real_,
               MAE = if (any(valid_mean)) mean(abs(d$observed_excess[valid_mean] - d$predicted_mean_excess[valid_mean])) else NA_real_,
               mean_conditional_log_score = mean(d$log_score), stringsAsFactors = FALSE)
  }))
}

# ==============================================================================
# STEP 3: SHARED PROJECT CONFIGURATION
# Source file: R/project_config.R
# ==============================================================================
# Single source of truth for the analysis and its knitted report.

project_config <- function() {
  list(
    threshold = 2,
    rolling_observations = 1250L,
    minimum_training_exceedances = 40L,
    forecast_step = 1L,
    seed = 20261018L,
    prior_specification = gpd_prior_specification()
  )
}

# ==============================================================================
# STEP 4: INDEPENDENT RUNTIME VALIDATION FUNCTIONS
# Source file: R/validation_functions.R
# ==============================================================================
# Independent checks for the data, timing, likelihood, rolling windows and outputs.

validate_analysis <- function(merged, predictions, diagnostics, metrics, input_path, cfg) {
  checks <- list()
  add_check <- function(name, passed, detail) {
    checks[[length(checks) + 1L]] <<- data.frame(
      check = name, passed = isTRUE(passed), detail = as.character(detail),
      stringsAsFactors = FALSE)
  }
  predictors <- c("nifty_vol_lag1", "sp500_vol_lag1", "gold_vol_lag1",
                  "india_vix_lag1", "usdinr_vol_lag1")
  prior_spec <- if (is.null(cfg$prior_specification)) gpd_prior_specification() else cfg$prior_specification

  add_check("Nifty dates unique and sorted",
            !anyNA(merged$Date) && !anyDuplicated(merged$Date) &&
              identical(order(merged$Date), seq_len(nrow(merged))),
            sprintf("%d forecast dates", nrow(merged)))
  expected_event <- as.integer(is.finite(merged$nifty_return) & merged$nifty_return < -cfg$threshold)
  add_check("Threshold indicator and excess response",
            identical(merged$exceedance, expected_event) &&
              all(merged$excess[expected_event == 1L] > 0) &&
              all(is.na(merged$excess[expected_event == 0L])),
            sprintf("threshold %.3f percentage points; %d exceedances",
                    cfg$threshold, sum(expected_event)))

  # Test the strict-before rule on a holiday/weekend-style toy calendar.
  toy_target <- as.Date(c("2024-01-08", "2024-01-09", "2024-01-10"))
  toy_source <- as.Date(c("2024-01-05", "2024-01-08", "2024-01-12"))
  toy <- latest_prior(toy_target, toy_source, c(5, 8, 12))
  add_check("Strictly prior-date synthetic calendar",
            identical(unname(toy), c(5, 8, 8)),
            "A source observation on the target date is excluded; weekends use the prior available source date.")

  # Compare the production lag function against a separate brute-force lookup
  # on evenly spaced rows from each actual source calendar.
  source_sheets <- read_market_sheets(input_path)
  date_column <- function(sheet) if ("DATE" %in% names(sheet)) "DATE" else "Date"
  lag_map <- list(
    nifty_vol_lag1 = list(sheet = source_sheets$nifty, date = date_column(source_sheets$nifty), value = "Daily empirical volatility (%)"),
    sp500_vol_lag1 = list(sheet = source_sheets$sp500, date = date_column(source_sheets$sp500), value = "Daily empirical volatility (%)"),
    gold_vol_lag1 = list(sheet = source_sheets$gold, date = date_column(source_sheets$gold), value = "Daily empirical volatility (%)"),
    india_vix_lag1 = list(sheet = source_sheets$vix, date = date_column(source_sheets$vix), value = "India VIX level (points)"),
    usdinr_vol_lag1 = list(sheet = source_sheets$usdinr, date = date_column(source_sheets$usdinr), value = "Daily empirical volatility (%)")
  )
  sample_rows <- unique(as.integer(round(seq(1, nrow(merged), length.out = min(100L, nrow(merged))))))
  for (nm in names(lag_map)) {
    m <- lag_map[[nm]]
    source_dates <- as.Date(m$sheet[[m$date]])
    source_values <- as.numeric(m$sheet[[m$value]])
    comparisons <- vapply(sample_rows, function(i) {
      available <- which(!is.na(source_dates) & source_dates < merged$Date[i] & is.finite(source_values))
      if (!length(available)) return(is.na(merged[[nm]][i]))
      latest <- available[which.max(source_dates[available])]
      is.finite(merged[[nm]][i]) && abs(merged[[nm]][i] - source_values[latest]) < 1e-10
    }, logical(1))
    add_check(paste("Actual strictly prior alignment:", nm), all(comparisons),
              sprintf("independent brute-force comparison on %d evenly spaced forecast dates", length(sample_rows)))
  }

  # Verify the exponential limit, a normalized positive-shape density, and
  # rejection outside the finite support for negative shape.
  z <- c(0.1, 0.7, 1.8)
  sigma <- 1.4
  exp_limit <- gpd_log_density(z, rep(log(sigma), length(z)), 0)
  exact_exp <- stats::dexp(z, rate = 1 / sigma, log = TRUE)
  add_check("GPD exponential limit",
            isTRUE(all.equal(exp_limit, exact_exp, tolerance = 1e-12)),
            "xi = 0 agrees with the exponential density of scale sigma.")
  integral <- tryCatch(stats::integrate(function(q) {
    exp(gpd_log_density(q, rep(log(sigma), length(q)), 0.2))
  }, lower = 0, upper = Inf, rel.tol = 1e-8)$value, error = function(e) NA_real_)
  add_check("GPD positive-shape density integrates to one",
            is.finite(integral) && abs(integral - 1) < 1e-7,
            sprintf("numerical integral %.12f", integral))
  support_rejected <- is.infinite(gpd_log_density(2.1, log(1), -0.5)) &&
    gpd_log_density(2.1, log(1), -0.5) < 0
  add_check("GPD finite support enforced",
            support_rejected && is.finite(gpd_log_density(1, log(1), -0.5)),
            "For xi = -0.5 and sigma = 1, z = 1 is inside support and z = 2.1 is rejected.")

  zero_penalties <- vapply(c("cauchy", "lasso", "ridge", "g-prior"), function(pr) {
    prior_penalty(c(0, 0, 0), pr, cbind(1, matrix(c(-1, 0, 1, 1, 0, -1), ncol = 2)), prior_spec)
  }, numeric(1))
  prior_design <- cbind(Intercept = 1, matrix(c(-1, 0, 1, 1, 0, -1), ncol = 2))
  prior_beta <- c(0, 0.2, -0.1)
  positive_g <- prior_penalty(prior_beta, "g-prior", prior_design, prior_spec)
  g_value <- if (identical(prior_spec$g_prior_g, "n")) nrow(prior_design) else as.numeric(prior_spec$g_prior_g)
  expected_g <- as.numeric(crossprod(prior_beta[-1L],
                                     crossprod(prior_design[, -1L, drop = FALSE]) %*% prior_beta[-1L])) /
    (2 * g_value)
  add_check("Prior penalty definitions",
            all(abs(zero_penalties) < 1e-12) && is.finite(positive_g) && positive_g > 0 &&
              isTRUE(all.equal(positive_g, expected_g, tolerance = 1e-12)),
            "All slope priors are centered at zero; the g-prior uses X'X / (2g), g = n.")

  gradient_design <- cbind(Intercept = rep(1, 5),
                           x1 = c(-1, -0.5, 0, 0.5, 1),
                           x2 = c(0.5, -1, 0, 1, -0.5),
                           x3 = c(1, 0.2, -0.7, 0.3, -0.8))
  gradient_z <- c(0.2, 0.35, 0.7, 1.1, 0.5)
  gradient_par <- c(log(0.8), 0.2, -0.15, 0.1, 0.18)
  gradient_errors <- vapply(c("cauchy", "lasso", "ridge", "g-prior"), function(pr) {
    analytic <- gpd_map_gradient(gradient_par, gradient_z, gradient_design, pr, prior_spec)
    numeric <- vapply(seq_along(gradient_par), function(j) {
      h <- 1e-6
      plus <- minus <- gradient_par
      plus[j] <- plus[j] + h
      minus[j] <- minus[j] - h
      (gpd_map_objective(plus, gradient_z, gradient_design, pr, prior_spec) -
         gpd_map_objective(minus, gradient_z, gradient_design, pr, prior_spec)) / (2 * h)
    }, numeric(1))
    max(abs(analytic - numeric))
  }, numeric(1))
  add_check("Analytic MAP gradients agree with finite differences",
            all(is.finite(gradient_errors)) && max(gradient_errors) < 1e-4,
            sprintf("maximum absolute difference %.3g across four priors", max(gradient_errors)))

  has_predictions <- nrow(predictions) > 0L
  add_check("Rolling predictions produced", has_predictions,
            sprintf("%d model-date rows", nrow(predictions)))
  if (has_predictions) {
    row_idx <- match(predictions$Date, merged$Date)
    event_ok <- !anyNA(row_idx) && all(merged$exceedance[row_idx] == 1L) &&
      isTRUE(all.equal(predictions$observed_excess, merged$excess[row_idx], tolerance = 1e-12)) &&
      all(complete.cases(merged[row_idx, predictors, drop = FALSE]))
    add_check("Forecast dates and responses match eligible events", event_ok,
              "Every scored response matches its source row and has complete strictly lagged covariates.")
    scales_ok <- all(is.finite(predictions$predicted_scale) & predictions$predicted_scale > 0)
    add_check("Predicted GPD scales are positive and finite", scales_ok,
              sprintf("%d predictions", nrow(predictions)))
    regression_rows <- predictions[predictions$model != "constant_gpd", , drop = FALSE]
    split_regression <- split(regression_rows, regression_rows$model)
    regression_models <- c("cauchy", "lasso", "ridge", "g-prior")
    distinct_scales <- all(regression_models %in% names(split_regression))
    if (distinct_scales) {
      common_dates <- Reduce(intersect, lapply(split_regression[regression_models], `[[`, "Date"))
      scale_matrix <- sapply(split_regression[regression_models], function(d) {
        d$predicted_scale[match(common_dates, d$Date)]
      })
      distinct_scales <- length(common_dates) > 0L &&
        any(apply(scale_matrix, 1L, function(v) max(v) - min(v) > 1e-8))
    }
    add_check("Prior-specific scales are not collapsed to one fit", distinct_scales,
              "At least one held-out event has different scales across the prior-regularized regressions.")
    score_ok <- !anyNA(predictions$log_score) && !any(is.nan(predictions$log_score)) &&
      !any(predictions$log_score == Inf)
    add_check("Held-out log-score values are valid", score_ok,
              sprintf("%d predictions outside fitted support", sum(predictions$log_score == -Inf)))
  }

  if (nrow(diagnostics)) {
    unique_fit_keys <- !anyDuplicated(paste(diagnostics$Date, diagnostics$model, sep = "|"))
    add_check("One diagnostic row per date and model", unique_fit_keys,
              sprintf("%d attempted fits", nrow(diagnostics)))
    successful <- !is.na(diagnostics$convergence) & diagnostics$convergence == 0L
    converged_ok <- all(is.finite(diagnostics$objective[successful])) &&
      all(nzchar(diagnostics$message[successful]))
    add_check("Converged fits have finite objectives and status messages", converged_ok,
              sprintf("%d successful, %d failed or skipped", sum(successful), sum(!successful)))
    regression_models <- c("cauchy", "lasso", "ridge", "g-prior")
    slope_evidence <- vapply(regression_models, function(model) {
      rows <- diagnostics$model == model & successful
      any(rows & is.finite(diagnostics$max_abs_slope) & diagnostics$max_abs_slope > 1e-8)
    }, logical(1))
    add_check("Prior regressions estimate nonzero slopes", all(slope_evidence),
              paste(names(slope_evidence), as.logical(slope_evidence), collapse = "; "))
    windows_ok <- vapply(seq_len(nrow(diagnostics)), function(j) {
      date <- diagnostics$Date[j]
      earlier <- sort(unique(merged$Date[merged$Date < date]))
      keep <- tail(earlier, cfg$rolling_observations)
      count <- sum(merged$Date %in% keep & merged$exceedance == 1L &
                     is.finite(merged$excess) & merged$excess > 0 &
                     complete.cases(merged[, predictors, drop = FALSE]))
      all(keep < date) && count == diagnostics$n_train[j]
    }, logical(1))
    add_check("Rolling-window dates and event counts", all(windows_ok),
              sprintf("recomputed preceding-window counts for %d diagnostic rows", length(windows_ok)))
  } else {
    add_check("Rolling fit diagnostics present", FALSE, "No eligible rolling fits were recorded.")
  }

  if (nrow(diagnostics) && has_predictions) {
    first_date <- min(diagnostics$Date)
    earlier <- sort(unique(merged$Date[merged$Date < first_date]))
    keep <- tail(earlier, cfg$rolling_observations)
    train <- merged[merged$Date %in% keep & merged$exceedance == 1L &
                      is.finite(merged$excess) & merged$excess > 0 &
                      complete.cases(merged[, predictors, drop = FALSE]), , drop = FALSE]
    target <- merged[merged$Date == first_date & merged$exceedance == 1L &
                       complete.cases(merged[, predictors, drop = FALSE]), , drop = FALSE]
    scaling_ok <- FALSE
    if (nrow(train) >= 2L && nrow(target)) {
      probe <- tryCatch(fit_gpd_map(train$excess, train[, predictors, drop = FALSE],
                                    "ridge", prior_spec),
                        error = function(e) NULL)
      if (!is.null(probe)) {
        expected_mean <- colMeans(train[, predictors, drop = FALSE])
        expected_sd <- apply(train[, predictors, drop = FALSE], 2L, stats::sd)
        manual_x <- cbind(Intercept = 1,
                          scale(as.matrix(target[, predictors, drop = FALSE]),
                                center = probe$means, scale = probe$sds))
        manual_scale <- exp(drop(manual_x %*% probe$beta))
        scaling_ok <- isTRUE(all.equal(unname(probe$means), unname(expected_mean), tolerance = 1e-12)) &&
          isTRUE(all.equal(unname(probe$sds), unname(expected_sd), tolerance = 1e-12)) &&
          isTRUE(all.equal(unname(predict(probe, target[, predictors, drop = FALSE])),
                           unname(manual_scale), tolerance = 1e-12))
      }
    }
    add_check("Training-only standardization and prediction reuse",
              scaling_ok,
              "Training means/SDs are calculated from training exceedances and reused for the held-out date.")
  }

  if (nrow(metrics) && has_predictions) {
    consistent <- all(vapply(seq_len(nrow(metrics)), function(j) {
      d <- predictions[predictions$model == metrics$model[j], , drop = FALSE]
      valid <- is.finite(d$predicted_mean_excess)
      rmse <- if (any(valid)) sqrt(mean((d$observed_excess[valid] - d$predicted_mean_excess[valid])^2)) else NA_real_
      mae <- if (any(valid)) mean(abs(d$observed_excess[valid] - d$predicted_mean_excess[valid])) else NA_real_
      isTRUE(metrics$forecasts[j] == nrow(d)) &&
        isTRUE(all.equal(metrics$RMSE[j], rmse, tolerance = 1e-12)) &&
        isTRUE(all.equal(metrics$MAE[j], mae, tolerance = 1e-12)) &&
        isTRUE(all.equal(metrics$mean_conditional_log_score[j], mean(d$log_score), tolerance = 1e-12))
    }, logical(1)))
    add_check("Metric table agrees with prediction rows", consistent,
              sprintf("independently recomputed metrics for %d models", nrow(metrics)))
  } else {
    add_check("Metric table present", FALSE, "No model metrics were generated.")
  }

  result <- do.call(rbind, checks)
  if (any(!result$passed)) {
    failed <- result$check[!result$passed]
    stop("Validation failed: ", paste(failed, collapse = "; "))
  }
  result
}

# ==============================================================================
# STEP 5: PRESPECIFIED REDUCED-MODEL COMPARISON (SETS A-D)
# Source file: R/reduced_model_comparison.R
# ==============================================================================
# Prespecified reduced predictor sets evaluated on common rolling dates.
# Source after R/model_functions.R, then call:
# run_reduced_model_comparison(merged, cfg, predictions, diagnostics)

reduced_model_sets <- function() {
  list(
    full_5 = c("nifty_vol_lag1", "sp500_vol_lag1", "gold_vol_lag1",
               "india_vix_lag1", "usdinr_vol_lag1"),
    A = c("nifty_vol_lag1", "sp500_vol_lag1", "gold_vol_lag1", "usdinr_vol_lag1"),
    B = c("sp500_vol_lag1", "gold_vol_lag1", "india_vix_lag1", "usdinr_vol_lag1"),
    C = c("nifty_vol_lag1", "usdinr_vol_lag1"),
    D = c("india_vix_lag1", "usdinr_vol_lag1")
  )
}

reduced_design_diagnostics <- function(x) {
  x <- as.data.frame(x)
  if (ncol(x) < 2L) return(c(max_abs_pairwise_r = NA_real_, max_vif = 1))
  cor_x <- stats::cor(x)
  max_abs_r <- max(abs(cor_x[upper.tri(cor_x)]))
  vif <- vapply(names(x), function(response) {
    other <- setdiff(names(x), response)
    fit <- stats::lm(stats::reformulate(other, response = response), data = x)
    r_squared <- summary(fit)$r.squared
    if (!is.finite(r_squared) || r_squared >= 1) Inf else 1 / (1 - r_squared)
  }, numeric(1))
  c(max_abs_pairwise_r = max_abs_r, max_vif = max(vif))
}

run_reduced_model_comparison <- function(merged, cfg, baseline_predictions,
                                         baseline_diagnostics, out_dir = "results") {
  full_predictors <- reduced_model_sets()$full_5
  model_sets <- reduced_model_sets()
  priors <- c("cauchy", "lasso", "ridge", "g-prior")
  required <- c("Date", "exceedance", "excess", full_predictors)
  if (!all(required %in% names(merged))) stop("Merged data is missing comparison columns")
  if (is.null(baseline_predictions) || !nrow(baseline_predictions) ||
      is.null(baseline_diagnostics) || !nrow(baseline_diagnostics)) {
    stop("Pass the existing rolling predictions and fit diagnostics to fix the evaluation dates")
  }

  dates <- sort(unique(as.Date(baseline_predictions$Date)))
  merged$Date <- as.Date(merged$Date)
  if (anyNA(dates) || !length(dates)) stop("No valid rolling forecast dates supplied")
  complete_full <- stats::complete.cases(merged[, full_predictors, drop = FALSE])
  eligible_dates <- merged$Date[merged$exceedance == 1L & is.finite(merged$excess) &
                                  merged$excess > 0 & complete_full]
  if (!all(dates %in% eligible_dates)) stop("A requested forecast date is not an eligible complete exceedance")

  prediction_rows <- list()
  fit_rows <- list()
  window_rows <- list()
  append_prediction <- function(row) {
    prediction_rows[[length(prediction_rows) + 1L]] <<- row
  }
  append_fit <- function(row) {
    fit_rows[[length(fit_rows) + 1L]] <<- row
  }

  for (date_value in dates) {
    # R's for-loop drops the Date class, so restore it before comparisons/exports.
    date <- as.Date(date_value, origin = "1970-01-01")
    earlier_dates <- sort(unique(merged$Date[merged$Date < date]))
    window_dates <- tail(earlier_dates, cfg$rolling_observations)
    if (length(window_dates) != cfg$rolling_observations || max(window_dates) >= date) {
      stop("A comparison window does not contain the configured number of strictly earlier dates")
    }
    window <- merged[merged$Date %in% window_dates, , drop = FALSE]
    complete_window <- stats::complete.cases(window[, full_predictors, drop = FALSE])
    train <- window[window$exceedance == 1L & is.finite(window$excess) & window$excess > 0 &
                      complete_window, , drop = FALSE]
    target <- merged[merged$Date == date & merged$exceedance == 1L &
                       is.finite(merged$excess) & merged$excess > 0 & complete_full, , drop = FALSE]
    if (nrow(train) < cfg$minimum_training_exceedances || nrow(target) != 1L) {
      stop("A comparison fold failed the training-event or single-target requirement")
    }
    window_rows[[length(window_rows) + 1L]] <- data.frame(
      forecast_date = date, window_first_date = min(window_dates),
      window_last_date = max(window_dates), n_training_exceedances = nrow(train),
      n_window_dates = length(window_dates), stringsAsFactors = FALSE
    )

    baseline <- tryCatch(fit_constant_gpd(train$excess, cfg$prior_specification),
                         error = function(e) e)
    if (inherits(baseline, "error")) stop("Constant GPD failed on ", date, ": ", conditionMessage(baseline))
    scale_hat <- exp(baseline$log_sigma)
    mean_hat <- if (baseline$xi < 1) scale_hat / (1 - baseline$xi) else NA_real_
    log_score <- gpd_log_density(target$excess, baseline$log_sigma, baseline$xi)
    append_prediction(data.frame(
      Date = date, model_set = "constant_gpd", prior = "baseline", predictor_count = 0L,
      predictors = "(intercept only)", observed_excess = target$excess,
      predicted_scale = scale_hat, predicted_mean_excess = mean_hat,
      shape = baseline$xi, log_score = log_score, stringsAsFactors = FALSE
    ))
    append_fit(data.frame(
      forecast_date = date, model_set = "constant_gpd", prior = "baseline",
      n_training_exceedances = nrow(train), predictor_count = 0L,
      max_abs_pairwise_r = NA_real_, max_vif = NA_real_,
      convergence = baseline$convergence, objective = baseline$objective,
      iterations = baseline$iterations, max_abs_slope = 0, error = "",
      stringsAsFactors = FALSE
    ))

    for (set_name in names(model_sets)) {
      selected <- model_sets[[set_name]]
      x_train <- train[, selected, drop = FALSE]
      x_target <- target[, selected, drop = FALSE]
      design_stats <- reduced_design_diagnostics(x_train)
      for (prior_name in priors) {
        fitted <- tryCatch(fit_gpd_map(train$excess, x_train, prior_name,
                                      cfg$prior_specification), error = function(e) e)
        if (inherits(fitted, "error")) {
          append_fit(data.frame(
            forecast_date = date, model_set = set_name, prior = prior_name,
            n_training_exceedances = nrow(train), predictor_count = length(selected),
            max_abs_pairwise_r = unname(design_stats["max_abs_pairwise_r"]),
            max_vif = unname(design_stats["max_vif"]), convergence = NA_integer_,
            objective = NA_real_, iterations = NA_integer_, max_abs_slope = NA_real_,
            error = conditionMessage(fitted), stringsAsFactors = FALSE
          ))
          next
        }
        predicted_scale <- predict(fitted, x_target)
        predicted_mean <- if (fitted$xi < 1) predicted_scale / (1 - fitted$xi) else NA_real_
        score <- gpd_log_density(target$excess, log(predicted_scale), fitted$xi)
        append_prediction(data.frame(
          Date = date, model_set = set_name, prior = prior_name,
          predictor_count = length(selected), predictors = paste(selected, collapse = ";"),
          observed_excess = target$excess, predicted_scale = predicted_scale,
          predicted_mean_excess = predicted_mean, shape = fitted$xi,
          log_score = score, stringsAsFactors = FALSE
        ))
        append_fit(data.frame(
          forecast_date = date, model_set = set_name, prior = prior_name,
          n_training_exceedances = nrow(train), predictor_count = length(selected),
          max_abs_pairwise_r = unname(design_stats["max_abs_pairwise_r"]),
          max_vif = unname(design_stats["max_vif"]),
          convergence = fitted$convergence, objective = fitted$objective,
          iterations = fitted$iterations, max_abs_slope = max(abs(fitted$beta[-1L])),
          error = "", stringsAsFactors = FALSE
        ))
      }
    }
  }

  predictions <- do.call(rbind, prediction_rows)
  diagnostics <- do.call(rbind, fit_rows)
  windows <- do.call(rbind, window_rows)
  expected_specs <- c("constant_gpd", names(model_sets))
  expected_count <- length(dates) * (1L + length(model_sets) * length(priors))
  groups <- split(predictions, interaction(predictions$model_set, predictions$prior, drop = TRUE))
  metrics <- do.call(rbind, lapply(groups, function(d) {
    valid_mean <- is.finite(d$predicted_mean_excess)
    data.frame(
      model_set = d$model_set[1], prior = d$prior[1],
      predictor_count = d$predictor_count[1], predictors = d$predictors[1],
      forecasts = nrow(d), finite_mean_forecasts = sum(valid_mean),
      RMSE = if (any(valid_mean)) sqrt(mean((d$observed_excess[valid_mean] -
                                               d$predicted_mean_excess[valid_mean])^2)) else NA_real_,
      MAE = if (any(valid_mean)) mean(abs(d$observed_excess[valid_mean] -
                                           d$predicted_mean_excess[valid_mean])) else NA_real_,
      mean_conditional_log_score = mean(d$log_score),
      failed_fits = sum(diagnostics$model_set == d$model_set[1] &
                          diagnostics$prior == d$prior[1] & is.na(diagnostics$convergence)),
      stringsAsFactors = FALSE
    )
  }))
  metrics <- metrics[order(match(metrics$model_set, expected_specs),
                           match(metrics$prior, c("baseline", priors))), ]
  rownames(metrics) <- NULL

  rank_candidates <- function(values, prior_values) {
    ranks <- rep(NA_real_, length(values))
    for (prior_name in priors) {
      rows <- prior_values == prior_name
      ranks[rows] <- rank(values[rows], ties.method = "min", na.last = "keep")
    }
    ranks
  }
  for (metric_name in c("RMSE", "MAE")) {
    metrics[[paste0(metric_name, "_rank_within_prior")]] <- rank_candidates(
      metrics[[metric_name]], metrics$prior
    )
  }
  metrics$log_score_rank_within_prior <- rank_candidates(
    -metrics$mean_conditional_log_score, metrics$prior
  )

  # Validate the full-five fits and baseline against the already-run pipeline.
  common_equal <- function(x, y) {
    isTRUE(all.equal(unname(x), unname(y), tolerance = 1e-7, check.attributes = FALSE))
  }
  original <- baseline_predictions
  original$Date <- as.Date(original$Date)
  original_diagnostics <- baseline_diagnostics
  original_diagnostics$Date <- as.Date(original_diagnostics$Date)
  reproduction_checks <- list()
  for (set_name in c("constant_gpd", "full_5")) {
    for (prior_name in if (set_name == "constant_gpd") "baseline" else priors) {
      model_name <- if (set_name == "constant_gpd") "constant_gpd" else prior_name
      new_rows <- predictions[predictions$model_set == set_name & predictions$prior == prior_name, ]
      old_rows <- original[original$model == model_name, ]
      new_rows <- new_rows[order(new_rows$Date), ]
      old_rows <- old_rows[order(old_rows$Date), ]
      same <- nrow(new_rows) == nrow(old_rows) && nrow(new_rows) == length(dates)
      if (same) {
        same <- identical(new_rows$Date, old_rows$Date) &&
          common_equal(new_rows$observed_excess, old_rows$observed_excess) &&
          common_equal(new_rows$predicted_scale, old_rows$predicted_scale) &&
          common_equal(new_rows$predicted_mean_excess, old_rows$predicted_mean_excess) &&
          common_equal(new_rows$shape, old_rows$shape) && common_equal(new_rows$log_score, old_rows$log_score)
      }
      reproduction_checks[[length(reproduction_checks) + 1L]] <- data.frame(
        check = paste("Matches original rolling pipeline:", set_name, prior_name),
        passed = isTRUE(same), detail = sprintf("%d comparison forecasts; %d original forecasts", nrow(new_rows), nrow(old_rows)),
        stringsAsFactors = FALSE
      )
    }
  }
  date_counts <- vapply(expected_specs, function(nm) {
    length(unique(predictions$Date[predictions$model_set == nm]))
  }, integer(1))
  combinations <- expand.grid(
    model_set = expected_specs,
    prior = c("baseline", priors),
    stringsAsFactors = FALSE
  )
  combinations <- combinations[
    (combinations$model_set == "constant_gpd" & combinations$prior == "baseline") |
      (combinations$model_set != "constant_gpd" & combinations$prior != "baseline"), , drop = FALSE
  ]
  combination_counts <- vapply(seq_len(nrow(combinations)), function(i) {
    sum(predictions$model_set == combinations$model_set[i] &
          predictions$prior == combinations$prior[i])
  }, integer(1))
  baseline_n <- original_diagnostics[original_diagnostics$model == "constant_gpd", c("Date", "n_train"), drop = FALSE]
  baseline_n <- baseline_n[order(baseline_n$Date), , drop = FALSE]
  window_n <- windows[order(windows$forecast_date), c("forecast_date", "n_training_exceedances"), drop = FALSE]
  same_training_counts <- nrow(baseline_n) == length(dates) && nrow(window_n) == length(dates) &&
    identical(baseline_n$Date, window_n$forecast_date) &&
    identical(as.integer(baseline_n$n_train), as.integer(window_n$n_training_exceedances))
  checks <- do.call(rbind, c(reproduction_checks, list(
    data.frame(check = "All requested models score every common forecast date",
               passed = nrow(predictions) == expected_count && all(date_counts == length(dates)) &&
                 all(combination_counts == length(dates)),
               detail = sprintf("%d rows expected/observed; %d dates for every set/prior combination",
                                expected_count, min(combination_counts)), stringsAsFactors = FALSE),
    data.frame(check = "Training event counts match the original five-predictor rolling folds",
               passed = same_training_counts,
               detail = sprintf("%d folds compared against original diagnostics", nrow(windows)),
               stringsAsFactors = FALSE),
    data.frame(check = "All model fits converged and predictions are finite/in support",
               passed = all(!is.na(diagnostics$convergence) & diagnostics$convergence == 0L) &&
                 all(is.finite(predictions$predicted_scale) & predictions$predicted_scale > 0) &&
                 all(is.finite(predictions$log_score)),
               detail = sprintf("%d converged fits; %d failed; %d nonfinite conditional scores",
                                sum(!is.na(diagnostics$convergence) & diagnostics$convergence == 0L),
                                sum(is.na(diagnostics$convergence) | diagnostics$convergence != 0L),
                                sum(!is.finite(predictions$log_score))), stringsAsFactors = FALSE),
    data.frame(check = "Training windows are strictly prior and meet the minimum event count",
               passed = all(windows$window_last_date < windows$forecast_date) &&
                 all(windows$n_window_dates == cfg$rolling_observations) &&
                 all(windows$n_training_exceedances >= cfg$minimum_training_exceedances),
               detail = sprintf("%d folds; %d-%d training exceedances", nrow(windows),
                                min(windows$n_training_exceedances), max(windows$n_training_exceedances)),
               stringsAsFactors = FALSE)
  )))

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  csv_copy <- function(x, path) {
    if ("Date" %in% names(x)) x$Date <- format(as.Date(x$Date), "%Y-%m-%d")
    if ("forecast_date" %in% names(x)) x$forecast_date <- format(as.Date(x$forecast_date), "%Y-%m-%d")
    if ("window_first_date" %in% names(x)) x$window_first_date <- format(as.Date(x$window_first_date), "%Y-%m-%d")
    if ("window_last_date" %in% names(x)) x$window_last_date <- format(as.Date(x$window_last_date), "%Y-%m-%d")
    utils::write.csv(x, file.path(out_dir, path), row.names = FALSE, na = "")
  }
  csv_copy(predictions, "reduced_model_comparison_predictions.csv")
  csv_copy(metrics, "reduced_model_comparison_metrics.csv")
  csv_copy(diagnostics, "reduced_model_comparison_fit_diagnostics.csv")
  csv_copy(windows, "reduced_model_comparison_windows.csv")
  csv_copy(checks, "reduced_model_comparison_validation.csv")
  saveRDS(list(model_sets = model_sets, predictions = predictions, metrics = metrics,
               diagnostics = diagnostics, windows = windows, validation = checks),
          file.path(out_dir, "reduced_model_comparison.rds"))

  if (any(!checks$passed)) {
    stop("Reduced comparison validation failed: ", paste(checks$check[!checks$passed], collapse = "; "))
  }
  list(model_sets = model_sets, predictions = predictions, metrics = metrics,
       diagnostics = diagnostics, windows = windows, validation = checks)
}

# ==============================================================================
# STEP 6: INTERACTIVE FILE/FOLDER POP-UPS & MAIN ANALYSIS PIPELINE (analysis.R)
# Opens pop-up dialogs to choose:
#   1. The input Excel workbook (RT_ONLY_COVARIATES.xlsx)
#   2. The output folder where results/ and presentation_materials/ will be saved
# ==============================================================================

select_input_file <- function() {
  message("Pop-up 1: Please select the input Excel file (RT_ONLY_COVARIATES.xlsx)...")
  path <- NULL
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    path <- rstudioapi::selectFile(
      caption = "Select Input Workbook (RT_ONLY_COVARIATES.xlsx)",
      filter = "Excel Files (*.xlsx)"
    )
  } else if (exists("choose.files", where = asNamespace("utils"), mode = "function")) {
    path <- utils::choose.files(
      caption = "Select Input Workbook (RT_ONLY_COVARIATES.xlsx)",
      multi = FALSE
    )
  } else {
    path <- base::file.choose()
  }
  if (is.null(path) || !length(path) || is.na(path) || !nzchar(path)) {
    stop("Execution stopped: No input Excel file was selected.")
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

select_output_folder <- function(default_dir = getwd()) {
  message("Pop-up 2: Please select the folder where you want to store the output files...")
  out_dir <- NULL
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    out_dir <- rstudioapi::selectDirectory(
      caption = "Select Folder to Store Output Data (results & presentation_materials)",
      path = default_dir
    )
  } else if (exists("choose.dir", where = asNamespace("utils"), mode = "function")) {
    out_dir <- utils::choose.dir(
      default = default_dir,
      caption = "Select Folder to Store Output Data (results & presentation_materials)"
    )
  } else if (requireNamespace("tcltk", quietly = TRUE)) {
    out_dir <- tcltk::tk_choose.dir(
      default = default_dir,
      caption = "Select Folder to Store Output Data (results & presentation_materials)"
    )
  } else {
    message("Select any file inside your target output folder to set the output directory:")
    out_dir <- dirname(base::file.choose())
  }
  if (is.null(out_dir) || !length(out_dir) || is.na(out_dir) || !nzchar(out_dir)) {
    stop("Execution stopped: No output folder was selected.")
  }
  normalizePath(out_dir, winslash = "/", mustWork = TRUE)
}

# 1. Choose input Excel file and output storage directory.
# [FIX B-13] Non-interactive runs (Rscript, CI) can pass paths instead of pop-ups:
#   Rscript pipeline.R <input.xlsx> <output_dir>   or set NIFTY_INPUT / NIFTY_OUTPUT.
cli_args <- commandArgs(trailingOnly = TRUE)
env_in <- Sys.getenv("NIFTY_INPUT"); env_out <- Sys.getenv("NIFTY_OUTPUT")
input_file <- if (length(cli_args) >= 1L) normalizePath(cli_args[1], mustWork = TRUE) else
  if (nzchar(env_in)) normalizePath(env_in, mustWork = TRUE) else select_input_file()
output_base_dir <- if (length(cli_args) >= 2L) { dir.create(cli_args[2], FALSE, TRUE); normalizePath(cli_args[2]) } else
  if (nzchar(env_out)) { dir.create(env_out, FALSE, TRUE); normalizePath(env_out) } else
  select_output_folder(default_dir = dirname(input_file))

results_dir <- file.path(output_base_dir, "results")
presentation_dir <- file.path(output_base_dir, "presentation_materials")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(presentation_dir, showWarnings = FALSE, recursive = TRUE)

message("Input workbook selected: ", input_file)
message("Outputs will be stored in: ", output_base_dir)

# 2. Run main analysis pipeline
cfg <- project_config()
merged <- prepare_and_merge(input_file, threshold = cfg$threshold)
set.seed(cfg$seed)
predictions <- rolling_evaluation(merged, cfg)
metrics <- score_predictions(predictions)
diagnostics <- attr(predictions, "fit_diagnostics")
reduced_comparison <- run_reduced_model_comparison(
  merged, cfg, predictions, diagnostics, out_dir = results_dir
)

# Export dates as ISO strings so CSV files are portable and unambiguous.
write_project_csv <- function(x, path) {
  if ("Date" %in% names(x)) x$Date <- format(x$Date, "%Y-%m-%d")
  utils::write.csv(x, path, row.names = FALSE, na = "")
}
# [FIX B-03] keep untouched copies; later steps re-use the names merged/predictions/metrics.
main_merged <- merged; main_predictions <- predictions; main_metrics <- metrics
write_project_csv(merged, file.path(results_dir, "merged_data.csv"))
write_project_csv(build_extremes(merged), file.path(results_dir, "extreme_severity_data.csv"))
write_project_csv(predictions, file.path(results_dir, "rolling_predictions.csv"))
write_project_csv(metrics, file.path(results_dir, "rolling_metrics.csv"))
write_project_csv(diagnostics, file.path(results_dir, "fit_diagnostics.csv"))

validation <- tryCatch(
  validate_analysis(merged, predictions, diagnostics, metrics, input_file, cfg),
  error = function(e) e
)
if (inherits(validation, "error")) {
  writeLines(conditionMessage(validation), file.path(results_dir, "validation_failure.txt"))
  stop(validation)
}
if (file.exists(file.path(results_dir, "validation_failure.txt"))) {
  unlink(file.path(results_dir, "validation_failure.txt"))
}
write_project_csv(validation, file.path(results_dir, "validation_checks.csv"))
session_info <- utils::sessionInfo()
writeLines(capture.output(session_info), file.path(results_dir, "session_info.txt"))
saveRDS(list(config = cfg, predictions = predictions, metrics = metrics,
             diagnostics = diagnostics, validation = validation,
             reduced_model_comparison = reduced_comparison,
             session_info = session_info),
        file.path(results_dir, "analysis_objects.rds"))
message("Analysis and validation finished. Outputs are saved in: ", results_dir)
message(sprintf("%d exceedances; %d rolling predictions; %d/%d validation checks passed.",
                sum(merged$exceedance == 1L), nrow(predictions), sum(validation$passed), nrow(validation)))
message(sprintf("Reduced-model comparison: %d prediction rows; %d/%d validation checks passed.",
                nrow(reduced_comparison$predictions), sum(reduced_comparison$validation$passed),
                nrow(reduced_comparison$validation)))

# ==============================================================================
# STEP 7: PREDICTOR MULTICOLLINEARITY DIAGNOSTICS
# Source file: R/collinearity_diagnostics.R
# ==============================================================================
# Read-only diagnostics for predictor collinearity in the fitted rolling design.
# Run from the project root with source("R/collinearity_diagnostics.R").


cfg <- project_config()
predictors <- c("nifty_vol_lag1", "sp500_vol_lag1", "gold_vol_lag1",
                "india_vix_lag1", "usdinr_vol_lag1")
merged <- utils::read.csv(file.path(results_dir, "merged_data.csv"), stringsAsFactors = FALSE)
predictions <- utils::read.csv(file.path(results_dir, "rolling_predictions.csv"), stringsAsFactors = FALSE)
merged$Date <- as.Date(merged$Date)
predictions$Date <- as.Date(predictions$Date)

if (!all(c("Date", "exceedance", "excess", predictors) %in% names(merged))) {
  stop("Merged results are missing required collinearity fields")
}

profile_design <- function(data) {
  x <- as.data.frame(lapply(data[, predictors, drop = FALSE], as.numeric))
  if (nrow(x) < 2L || any(!vapply(x, function(z) is.finite(stats::sd(z)), logical(1))) ||
      any(vapply(x, function(z) stats::sd(z) == 0, logical(1)))) {
    stop("Collinearity design has too few rows or a zero-variance predictor")
  }

  correlation <- stats::cor(x, method = "pearson")
  upper <- which(upper.tri(correlation), arr.ind = TRUE)
  pair_values <- correlation[upper]
  max_pair_index <- which.max(abs(pair_values))
  max_pair <- paste(rownames(correlation)[upper[max_pair_index, 1]],
                    colnames(correlation)[upper[max_pair_index, 2]], sep = " -- ")

  vif <- vapply(seq_along(predictors), function(j) {
    response <- predictors[j]
    other <- predictors[-j]
    fit <- stats::lm(stats::reformulate(other, response = response), data = x)
    r_squared <- summary(fit)$r.squared
    if (!is.finite(r_squared) || r_squared >= 1) Inf else 1 / (1 - r_squared)
  }, numeric(1))
  names(vif) <- paste0("vif_", predictors)

  eigenvalues <- eigen(correlation, symmetric = TRUE, only.values = TRUE)$values
  smallest <- min(eigenvalues)
  condition_index <- if (smallest <= 0) Inf else sqrt(max(eigenvalues) / smallest)

  list(
    summary = data.frame(
      training_exceedances = nrow(x),
      design_rank = qr(cbind(`(Intercept)` = 1, as.matrix(x)))$rank,
      max_abs_pairwise_r = abs(pair_values[max_pair_index]),
      max_correlated_pair = max_pair,
      max_vif = max(vif),
      max_vif_predictor = sub("^vif_", "", names(vif)[which.max(vif)]),
      condition_index = condition_index,
      stringsAsFactors = FALSE
    ),
    correlation = correlation,
    vif = vif
  )
}

# Full-sample diagnostic is descriptive. The rolling diagnostics below mirror
# each model fit: preceding Nifty dates, then complete training exceedances.
full_events <- merged[
  merged$exceedance == 1L & is.finite(merged$excess) & merged$excess > 0 &
    stats::complete.cases(merged[, predictors, drop = FALSE]), , drop = FALSE
]
full_profile <- profile_design(full_events)
utils::write.csv(full_profile$correlation,
                 file.path(results_dir, "collinearity_full_sample_pearson.csv"), row.names = TRUE)
utils::write.csv(stats::cor(full_events[, predictors, drop = FALSE], method = "spearman"),
                 file.path(results_dir, "collinearity_full_sample_spearman.csv"), row.names = TRUE)
full_vif <- data.frame(
  predictor = predictors,
  VIF = as.numeric(full_profile$vif),
  stringsAsFactors = FALSE
)
utils::write.csv(full_vif, file.path(results_dir, "collinearity_full_sample_vif.csv"), row.names = FALSE)

forecast_dates <- sort(unique(predictions$Date))
if (!length(forecast_dates)) stop("No held-out event dates found in rolling predictions")

rolling_rows <- vector("list", length(forecast_dates))
pair_rows <- vector("list", length(forecast_dates))
for (i in seq_along(forecast_dates)) {
  target_date <- forecast_dates[i]
  prior_dates <- sort(unique(merged$Date[merged$Date < target_date]))
  window_dates <- tail(prior_dates, cfg$rolling_observations)
  if (length(window_dates) != cfg$rolling_observations ||
      max(window_dates) >= target_date) {
    stop("Rolling diagnostic window does not match the configured prior-date window")
  }
  window_data <- merged[merged$Date %in% window_dates, , drop = FALSE]
  training_events <- window_data[
    window_data$exceedance == 1L & is.finite(window_data$excess) & window_data$excess > 0 &
      stats::complete.cases(window_data[, predictors, drop = FALSE]), , drop = FALSE
  ]
  if (nrow(training_events) < cfg$minimum_training_exceedances) {
    stop("Rolling collinearity window has fewer than the configured minimum events")
  }

  p <- profile_design(training_events)
  rolling_rows[[i]] <- cbind(
    data.frame(forecast_date = target_date,
               window_first_date = min(window_dates),
               window_last_date = max(window_dates),
               stringsAsFactors = FALSE),
    p$summary,
    as.data.frame(as.list(p$vif), check.names = FALSE)
  )

  upper <- which(upper.tri(p$correlation), arr.ind = TRUE)
  pair_rows[[i]] <- data.frame(
    forecast_date = target_date,
    predictor_1 = rownames(p$correlation)[upper[, 1]],
    predictor_2 = colnames(p$correlation)[upper[, 2]],
    pearson_r = p$correlation[upper],
    stringsAsFactors = FALSE
  )
}

rolling <- do.call(rbind, rolling_rows)
pairs <- do.call(rbind, pair_rows)
if (any(rolling$design_rank != length(predictors) + 1L) ||
    any(rolling$training_exceedances < cfg$minimum_training_exceedances) ||
    any(rolling$window_last_date >= rolling$forecast_date)) {
  stop("A rolling collinearity check failed rank, sample-size, or date-order validation")
}
utils::write.csv(rolling, file.path(results_dir, "collinearity_rolling_windows.csv"), row.names = FALSE)
utils::write.csv(pairs, file.path(results_dir, "collinearity_rolling_pairwise.csv"), row.names = FALSE)

pair_summary <- do.call(rbind, lapply(split(pairs, interaction(pairs$predictor_1, pairs$predictor_2,
                                                               drop = TRUE)), function(d) {
  data.frame(
    predictor_1 = d$predictor_1[1],
    predictor_2 = d$predictor_2[1],
    full_sample_pearson_r = full_profile$correlation[d$predictor_1[1], d$predictor_2[1]],
    rolling_median_r = stats::median(d$pearson_r),
    rolling_min_r = min(d$pearson_r),
    rolling_max_r = max(d$pearson_r),
    rolling_max_abs_r = max(abs(d$pearson_r)),
    stringsAsFactors = FALSE
  )
}))
pair_summary <- pair_summary[order(-pair_summary$rolling_max_abs_r), ]
utils::write.csv(pair_summary, file.path(results_dir, "collinearity_pair_summary.csv"), row.names = FALSE)

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
utils::write.csv(summary_row, file.path(results_dir, "collinearity_summary.csv"), row.names = FALSE)

# Produce a concise human-readable diagnostic alongside the machine-readable
# tables. Thresholds below are descriptive rules of thumb, not formal tests.
rolling_vif_summary <- do.call(rbind, lapply(predictors, function(variable) {
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
full_pair_summary <- data.frame(
  predictor_1 = rownames(full_profile$correlation)[full_upper[, 1]],
  predictor_2 = colnames(full_profile$correlation)[full_upper[, 2]],
  pearson_r = full_profile$correlation[full_upper],
  stringsAsFactors = FALSE
)
full_pair_summary <- full_pair_summary[order(-abs(full_pair_summary$pearson_r)), ]
top_full_pairs <- utils::head(full_pair_summary, 5)
top_rolling_pairs <- utils::head(pair_summary, 5)
vif_lines <- c(
  "| Predictor | Full sample VIF | Rolling median VIF | Rolling maximum VIF | Folds VIF >= 5 | Folds VIF >= 10 |",
  "|:--|--:|--:|--:|--:|--:|",
  vapply(seq_len(nrow(rolling_vif_summary)), function(i) {
    variable <- rolling_vif_summary$predictor[i]
    full_value <- full_vif$VIF[match(variable, full_vif$predictor)]
    sprintf("| %s | %.2f | %.2f | %.2f | %d/%d | %d/%d |",
            variable, full_value, rolling_vif_summary$median_vif[i],
            rolling_vif_summary$maximum_vif[i],
            rolling_vif_summary$windows_vif_ge_5[i], nrow(rolling),
            rolling_vif_summary$windows_vif_ge_10[i], nrow(rolling))
  }, character(1))
)
full_pair_lines <- c(
  "| Predictor pair | Pearson r |",
  "|:--|--:|",
  vapply(seq_len(nrow(top_full_pairs)), function(i) {
    sprintf("| %s / %s | %.3f |", top_full_pairs$predictor_1[i],
            top_full_pairs$predictor_2[i], top_full_pairs$pearson_r[i])
  }, character(1))
)
rolling_pair_lines <- c(
  "| Predictor pair | Median rolling r | Maximum rolling absolute r |",
  "|:--|--:|--:|",
  vapply(seq_len(nrow(top_rolling_pairs)), function(i) {
    sprintf("| %s / %s | %.3f | %.3f |", top_rolling_pairs$predictor_1[i],
            top_rolling_pairs$predictor_2[i], top_rolling_pairs$rolling_median_r[i],
            top_rolling_pairs$rolling_max_abs_r[i])
  }, character(1))
)
report_lines <- c(
  "# Predictor multicollinearity check",
  "",
  "This report audits the five lagged predictors on complete Nifty exceedance observations, matching the severity model's data. The rolling audit reconstructs the model's preceding-date windows and checks every held-out forecast date.",
  "",
  sprintf("The full-sample diagnostic used %d complete exceedances. The rolling audit covered %d forecast dates, with %d to %d training exceedances per fold and %d preceding Nifty dates per window.",
          nrow(full_events), nrow(rolling), min(rolling$training_exceedances),
          max(rolling$training_exceedances), cfg$rolling_observations),
  "",
  "## Main findings",
  "",
  sprintf("- The strongest full-sample Pearson association was **%s / %s (r = %.3f)**.",
          top_full_pairs$predictor_1[1], top_full_pairs$predictor_2[1], top_full_pairs$pearson_r[1]),
  sprintf("- Across rolling folds, the median of each fold's strongest absolute pairwise correlation was **%.3f**, and the largest was **%.3f**; %d of %d folds had a strongest pair with |r| >= 0.8.",
          summary_row$median_rolling_max_abs_r, summary_row$maximum_rolling_max_abs_r,
          summary_row$windows_max_abs_r_ge_0_8, nrow(rolling)),
  sprintf("- The median fold maximum VIF was **%.2f**, the maximum was **%.2f**; %d of %d folds had a maximum VIF >= 5 and %d had a maximum VIF >= 10.",
          summary_row$median_rolling_max_vif, summary_row$maximum_rolling_vif,
          summary_row$windows_max_vif_ge_5, nrow(rolling),
          summary_row$windows_max_vif_ge_10),
  sprintf("- All rolling model designs had rank %d (intercept plus all five predictors estimable); the maximum scaled-correlation condition index was **%.2f**.",
          length(predictors) + 1L, summary_row$maximum_rolling_condition_index),
  "",
  "A VIF around 5 is commonly treated as a warning and 10 as a stronger warning; these are heuristics, not pass/fail tests. Correlation and VIF are invariant to linear rescaling, so calculating them on the raw lagged predictors gives the same collinearity result as training-only standardisation.",
  "",
  "## Full-sample VIF by predictor",
  "",
  vif_lines,
  "",
  "## Largest full-sample pairwise correlations",
  "",
  full_pair_lines,
  "",
  "## Largest associations seen in rolling folds",
  "",
  rolling_pair_lines,
  "",
  "## Is this evidence of overfitting?",
  "",
  "It is evidence of redundant predictors and coefficient-instability risk, especially for Nifty volatility, S&P 500 volatility, and India VIX. It does not by itself establish that forecasts are overfit. At the smallest fold, 40 exceedances support five scale slopes, an intercept, and a shape estimate, so the effective information per fitted quantity is limited even though priors regularise the slopes.",
  "",
  # [FIX B-05] numbers computed from the saved metrics instead of hard-coded
  with(list(rm = utils::read.csv(file.path(results_dir, "rolling_metrics.csv"))), {
    reg <- rm[rm$model != "constant_gpd", ]; cst <- rm[rm$model == "constant_gpd", ]
    sprintf("The existing %d-date rolling results are mixed: covariate models have severity RMSE %.3f-%.3f versus %.3f for the constant GPD baseline, while their mean conditional log scores are %.3f to %.3f versus %.3f (higher is better). This is a reason to scrutinise predictive distribution calibration, but %d forecast events are too few to make a strong generalisation claim. The scores are conditional on an exceedance and do not assess the probability that an exceedance occurs.",
            nrow(rolling), min(reg$RMSE), max(reg$RMSE), cst$RMSE, max(reg$mean_conditional_log_score),
            min(reg$mean_conditional_log_score), cst$mean_conditional_log_score, nrow(rolling))
  }),
  "",
  "## Recommended follow-up",
  "",
  "1. Compare the current five-predictor fit with prespecified reduced sets, especially Nifty realised volatility versus India VIX and S&P 500 volatility versus gold volatility. Do not remove variables based on correlation alone; decide which signal is more useful and available at forecast time.",
  "2. Refit each reduced set on the same rolling folds and compare RMSE, MAE, and conditional log score. If selecting predictors or tuning prior strength, do that inside a nested rolling validation loop, keeping the final evaluation dates untouched.",
  "3. Treat the ridge and lasso fits as regularised sensitivity checks; consider increasing the minimum number of training exceedances or using a wider rolling window if that matches the intended regime assumptions.",
  "",
  "## Reproduce",
  "",
  "From the project root, run `source(\"R/collinearity_diagnostics.R\")` in R. Outputs are written to `results/collinearity_*.csv` and this report is rewritten as `results/collinearity_report.md`."
)
writeLines(report_lines, file.path(results_dir, "collinearity_report.md"), useBytes = TRUE)

cat("Full-sample complete exceedances:", nrow(full_events), "\n")
cat("Rolling windows audited:", nrow(rolling), "\n")
cat("Full-sample VIFs:\n")
print(full_vif, row.names = FALSE)
cat("Rolling summary:\n")
print(summary_row, row.names = FALSE)
cat("Largest rolling pairwise correlations:\n")
print(utils::head(pair_summary, 5), row.names = FALSE)

# ==============================================================================
# STEP 8: PRESENTATION TABLES & CHARTS GENERATOR
# Source file: R/create_presentation_materials.R
# ==============================================================================
# Create presentation-ready charts and tables from the saved project results.
# Run from the project root: source("R/create_presentation_materials.R")

suppressPackageStartupMessages(library(ggplot2))

results_dir <- file.path(output_base_dir, "results")
presentation_dir <- file.path(output_base_dir, "presentation_materials")
if (!dir.exists(presentation_dir)) dir.create(presentation_dir, recursive = TRUE)

read_project_csv <- function(file) {
  path <- file.path(results_dir, file)
  if (!file.exists(path)) stop("Required project result is missing: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

nice_names <- c(
  constant_gpd = "Constant GPD baseline (no covariates)",
  full_5 = "Nifty volatility + S&P 500 volatility + Gold volatility + India VIX + USD/INR volatility",
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

theme_presentation <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.2)),
      plot.subtitle = element_text(color = "#4B5563"),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title = element_blank()
    )
}

save_plot <- function(plot, filename, width = 11, height = 6) {
  ggsave(
    filename = file.path(presentation_dir, filename),
    plot = plot, width = width, height = height, units = "in",
    dpi = 300, bg = "white"
  )
}

# Table 1: data and evaluation design.
merged <- read_project_csv("merged_data.csv")
merged$Date <- as.Date(merged$Date)
event_count <- sum(as.numeric(merged$exceedance) == 1, na.rm = TRUE)
forecast_predictions <- read_project_csv("reduced_model_comparison_predictions.csv")
forecast_dates <- sort(unique(as.Date(forecast_predictions$Date)))
thresholds <- unique(na.omit(as.numeric(merged$threshold)))
data_summary <- data.frame(
  Item = c(
    "Nifty 50 data start", "Nifty 50 data end", "Extreme-loss threshold (percent)",
    "Extreme-loss exceedances", "Rolling forecast dates", "Prior Nifty trading dates per rolling window",
    "Forecast target", "Forecast scope"
  ),
  Value = c(
    format(min(merged$Date, na.rm = TRUE)), format(max(merged$Date, na.rm = TRUE)),
    paste(format(thresholds, trim = TRUE), collapse = ", "), event_count, length(forecast_dates),
    1250, "Excess loss above the threshold", "Severity conditional on threshold exceedance"
  ),
  stringsAsFactors = FALSE
)
write.csv(data_summary, file.path(presentation_dir, "table_data_and_evaluation_summary.csv"), row.names = FALSE)

# Table 2: all predictor-set/prior combinations with a clearly defined forecast R-squared.
metrics <- read_project_csv("reduced_model_comparison_metrics.csv")
fit_diagnostics <- read_project_csv("reduced_model_comparison_fit_diagnostics.csv")
prediction_groups <- split(
  forecast_predictions,
  interaction(forecast_predictions$model_set, forecast_predictions$prior, drop = TRUE)
)
forecast_r2_values <- vapply(prediction_groups, function(d) {
  actual <- as.numeric(d$observed_excess)
  predicted <- as.numeric(d$predicted_mean_excess)
  1 - sum((actual - predicted)^2) / sum((actual - mean(actual))^2)
}, numeric(1))
forecast_r2 <- data.frame(
  model_set_prior = names(forecast_r2_values),
  forecast_R_squared = as.numeric(forecast_r2_values),
  stringsAsFactors = FALSE
)
forecast_r2$model_set <- sub("\\..*$", "", forecast_r2$model_set_prior)
forecast_r2$prior <- sub("^.*\\.", "", forecast_r2$model_set_prior)

vif_rows <- unique(fit_diagnostics[c("forecast_date", "model_set", "max_vif")])
vif_rows$max_vif <- as.numeric(vif_rows$max_vif)
vif_rows <- vif_rows[vif_rows$model_set != "constant_gpd" & is.finite(vif_rows$max_vif), ]
vif_summary <- do.call(rbind, lapply(split(vif_rows, vif_rows$model_set), function(d) {
  data.frame(
    model_set = d$model_set[1],
    median_rolling_max_VIF = median(d$max_vif),
    maximum_rolling_VIF = max(d$max_vif),
    stringsAsFactors = FALSE
  )
}))

performance <- merge(metrics, forecast_r2[c("model_set", "prior", "forecast_R_squared")],
                      by = c("model_set", "prior"), all.x = TRUE, sort = FALSE)
performance <- merge(performance, vif_summary, by = "model_set", all.x = TRUE, sort = FALSE)
performance$Model <- unname(nice_names[performance$model_set])
performance$Prior <- performance$prior
performance$Forecasts <- performance$forecasts
presentation_metrics <- performance[c(
  "Model", "Prior", "Forecasts", "RMSE", "MAE", "forecast_R_squared",
  "mean_conditional_log_score", "median_rolling_max_VIF", "maximum_rolling_VIF"
)]
names(presentation_metrics) <- c(
  "Model", "Prior", "Forecast dates", "RMSE", "MAE", "Out-of-sample R-squared",
  "Mean conditional log score", "Median rolling max VIF", "Maximum rolling VIF"
)
write.csv(presentation_metrics, file.path(presentation_dir, "table_all_model_results.csv"), row.names = FALSE)

# Table 3: model definitions, spelled out for presentation use.
model_definitions <- data.frame(
  Model = unname(nice_names[c("constant_gpd", "full_5", "A", "B", "C", "D")]),
  Covariates = c(
    "None; intercept only",
    "Nifty volatility; S&P 500 volatility; Gold volatility; India VIX; USD/INR volatility",
    "Nifty volatility; S&P 500 volatility; Gold volatility; USD/INR volatility",
    "S&P 500 volatility; Gold volatility; India VIX; USD/INR volatility",
    "Nifty 50 volatility; USD/INR volatility",
    "India VIX; USD/INR volatility"
  ),
  stringsAsFactors = FALSE
)
write.csv(model_definitions, file.path(presentation_dir, "table_model_definitions.csv"), row.names = FALSE)

# Table 4: rolling collinearity summary by predictor set.
vif_summary$Model <- unname(nice_names[vif_summary$model_set])
vif_table <- vif_summary[c("Model", "median_rolling_max_VIF", "maximum_rolling_VIF")]
names(vif_table) <- c("Model", "Median rolling max VIF", "Maximum rolling VIF")
write.csv(vif_table, file.path(presentation_dir, "table_rolling_VIF_summary.csv"), row.names = FALSE)

# Graph 1: daily returns and threshold events.
return_data <- merged[!is.na(merged$nifty_return), c("Date", "nifty_return", "exceedance")]
return_data$exceedance <- as.numeric(return_data$exceedance)
threshold_value <- if (length(thresholds) == 1) thresholds else 2
p_returns <- ggplot(return_data, aes(Date, nifty_return)) +
  geom_hline(yintercept = -threshold_value, color = "#C0392B", linetype = "dashed", linewidth = 0.8) +
  geom_line(color = "#526D82", linewidth = 0.25, alpha = 0.8) +
  geom_point(data = return_data[return_data$exceedance == 1, ], color = "#C0392B", size = 1.1, alpha = 0.75) +
  labs(
    title = "Nifty 50 daily returns and extreme-loss threshold",
    subtitle = paste0(event_count, " threshold exceedances; red dashed line marks a loss of ", threshold_value, "%"),
    x = NULL, y = "Daily return (%)"
  ) + theme_presentation()
save_plot(p_returns, "graph_1_nifty_returns_threshold.png", 12, 5.5)

# Graph 2: observed and predicted severity on the same 28 rolling forecast dates.
chosen <- forecast_predictions[
  (forecast_predictions$model_set == "constant_gpd" & forecast_predictions$prior == "baseline") |
    (forecast_predictions$model_set == "C" & forecast_predictions$prior == "g-prior") |
    (forecast_predictions$model_set == "D" & forecast_predictions$prior == "cauchy"),
]
chosen$Date <- as.Date(chosen$Date)
chosen <- chosen[order(chosen$Date, chosen$model_set), ]
actual_line <- unique(chosen[c("Date", "observed_excess")])
actual_line$Series <- "Observed excess loss"
actual_line$Value <- actual_line$observed_excess
predicted_lines <- chosen[c("Date", "model_set", "prior", "predicted_mean_excess")]
predicted_lines$Series <- ifelse(
  predicted_lines$model_set == "constant_gpd", "Constant GPD baseline",
  ifelse(predicted_lines$model_set == "C",
         "Nifty volatility + USD/INR volatility (g-prior)",
         "India VIX + USD/INR volatility (Cauchy)")
)
predicted_lines$Value <- predicted_lines$predicted_mean_excess
plot_lines <- rbind(actual_line[c("Date", "Series", "Value")], predicted_lines[c("Date", "Series", "Value")])
p_forecast <- ggplot(plot_lines, aes(Date, Value, color = Series, group = Series)) +
  geom_line(linewidth = 0.75, alpha = 0.9) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c(
    "Observed excess loss" = "#111827", "Constant GPD baseline" = "#9CA3AF",
    "Nifty volatility + USD/INR volatility (g-prior)" = "#2563EB",
    "India VIX + USD/INR volatility (Cauchy)" = "#D97706"
  )) +
  labs(
    title = "Observed and predicted extreme-loss severity",
    subtitle = "Predictions are conditional mean excess losses for the same 28 rolling forecast dates",
    x = NULL, y = "Excess loss above threshold (%)"
  ) + theme_presentation(11) +
  theme(legend.text = element_text(size = 9))
save_plot(p_forecast, "graph_2_actual_vs_predicted_severity.png", 12, 6)

# Graph 3: fold-by-fold maximum VIF for each predictor set.
vif_rows$Model <- unname(nice_names[vif_rows$model_set])
vif_rows$Model <- factor(
  vif_rows$Model,
  levels = vif_table$Model[order(vif_table[["Maximum rolling VIF"]])]
)
p_vif <- ggplot(vif_rows, aes(max_vif, Model)) +
  geom_vline(xintercept = 5, linetype = "dashed", color = "#D97706", linewidth = 0.7) +
  geom_vline(xintercept = 10, linetype = "dotted", color = "#C0392B", linewidth = 0.7) +
  geom_boxplot(fill = "#8FB9A8", color = "#355C50", outlier.color = "#C0392B", outlier.size = 1.1) +
  labs(
    title = "Rolling multicollinearity by predictor set",
    subtitle = "Each box summarizes fold-specific maximum VIF; reference lines at 5 and 10 are common heuristics",
    x = "Maximum VIF within training fold", y = NULL
  ) + theme_presentation(10)
save_plot(p_vif, "graph_3_rolling_VIF_by_model.png", 12, 6.5)

# Graph 4: full-sample predictor correlation heatmap.
# [FIX B-02] use the 178 complete exceedance rows (the regression data), not all
# 4,390 complete trading days, so the heatmap matches the text and the VIF audit.
complete_x <- merged[as.numeric(merged$exceedance) == 1 & is.finite(merged$excess) & merged$excess > 0 &
                       complete.cases(merged[predictor_names]), predictor_names]
correlation <- cor(complete_x)
write.csv(correlation, file.path(presentation_dir, "table_predictor_correlation_matrix.csv"), row.names = TRUE)
cor_long <- as.data.frame(as.table(correlation), stringsAsFactors = FALSE)
names(cor_long) <- c("Predictor_1", "Predictor_2", "Correlation")
cor_long$Predictor_1 <- unname(predictor_labels[as.character(cor_long$Predictor_1)])
cor_long$Predictor_2 <- unname(predictor_labels[as.character(cor_long$Predictor_2)])
p_cor <- ggplot(cor_long, aes(Predictor_1, Predictor_2, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", Correlation)), size = 3.5, color = "#111827") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1, 1)) +
  coord_equal() +
  labs(
    title = "Correlation among lagged predictors",
    subtitle = paste0("Pearson correlations on the ", nrow(complete_x), " exceedances complete for all five covariates"),
    x = NULL, y = NULL, fill = "Correlation"
  ) + theme_presentation(10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_plot(p_cor, "graph_4_predictor_correlation_heatmap.png", 9, 7)

# Graph 5: compact error comparison for the baseline and the leading reduced specifications.
# The complete 21-row comparison is in table_all_model_results.csv.
shortlist <- performance[
  (performance$model_set == "constant_gpd" & performance$prior == "baseline") |
    (performance$model_set == "full_5" & performance$prior == "g-prior") |
    (performance$model_set == "C" & performance$prior == "g-prior") |
    (performance$model_set == "D" & performance$prior == "cauchy"),
]
shortlist$Specification <- ifelse(
  shortlist$model_set == "constant_gpd", "Constant GPD baseline",
  ifelse(shortlist$model_set == "full_5", "Full five-covariate model (g-prior)",
         ifelse(shortlist$model_set == "C", "Nifty volatility + USD/INR volatility (g-prior)",
                "India VIX + USD/INR volatility (Cauchy)"))
)
err_long <- rbind(
  data.frame(Specification = shortlist$Specification, Metric = "RMSE", Value = as.numeric(shortlist$RMSE)),
  data.frame(Specification = shortlist$Specification, Metric = "MAE", Value = as.numeric(shortlist$MAE))
)
err_long$Specification <- factor(err_long$Specification, levels = rev(shortlist$Specification))
p_errors <- ggplot(err_long, aes(Value, Specification, fill = Specification)) +
  geom_col(show.legend = FALSE, width = 0.68) +
  geom_text(aes(label = sprintf("%.3f", Value)), hjust = -0.15, size = 3.4) +
  facet_wrap(~Metric, scales = "free_x") +
  # [FIX B-12] name the colours so each specification keeps the colour used in Graph 2
  scale_fill_manual(values = c("Constant GPD baseline" = "#9CA3AF", "Full five-covariate model (g-prior)" = "#566573",
                               "Nifty volatility + USD/INR volatility (g-prior)" = "#2563EB",
                               "India VIX + USD/INR volatility (Cauchy)" = "#D97706")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(
    title = "Rolling forecast error for selected specifications",
    subtitle = "All values use the same 28 forecast dates; lower values indicate smaller errors",
    x = "Forecast error", y = NULL
  ) + theme_presentation(10) +
  theme(axis.text.y = element_text(size = 8))
save_plot(p_errors, "graph_5_selected_model_errors.png", 13, 6.5)

# Reproducibility and presentation notes.
notes <- c(
  "# Presentation materials",
  "",
  "Generated from the project's saved rolling forecasts and diagnostics by `R/create_presentation_materials.R`.",
  "",
  "## Files",
  "- Graph 1: daily Nifty 50 returns, threshold, and exceedance dates.",
  "- Graph 2: actual excess severity versus conditional mean forecasts for the baseline, Nifty volatility + USD/INR volatility (g-prior), and India VIX + USD/INR volatility (Cauchy).",
  "- Graph 3: distribution of fold-specific maximum VIF by covariate set; VIF reference values are heuristics.",
  "- Graph 4: full-sample Pearson correlation heatmap among the five lagged predictors.",
  "- Graph 5: RMSE and MAE comparison for selected specifications. The full prior-by-model results are in the results table.",
  "- Tables: event/evaluation summary, all model-prior results, model definitions, rolling VIF summary, and predictor correlation matrix.",
  "",
  "## Metric definitions and cautions",
  "- RMSE and MAE score predicted conditional mean excess loss against observed excess loss on the 28 forecast dates.",
  "- Out-of-sample R-squared is calculated as 1 - SSE/SST, where SST uses the mean observed excess over the 28 forecast dates. It is not adjusted R-squared and can be negative.",
  "- Conditional log score is averaged over forecast dates; a larger (less negative) value is better.",
  "- Rolling VIF is computed from each training fold. It diagnoses predictor collinearity; it does not itself establish predictive overfitting.",
  "- The GPD model predicts excess-loss severity conditional on a threshold exceedance, not the probability that an exceedance occurs.",
  "- The 28 dates were used for model comparison and selection, so the displayed winner should be treated as exploratory pending a separate untouched evaluation."
)
writeLines(notes, file.path(presentation_dir, "README.md"))

message("Presentation materials created in: ", normalizePath(presentation_dir, winslash = "/"))
message("Forecast dates: ", length(forecast_dates), "; exceedances: ", event_count,
        "; model/prior result rows: ", nrow(presentation_metrics))

# ==============================================================================
# STEP 9: REPORT TABLES & DIAGNOSTIC PLOTS (from report.Rmd)
# Generates the summary tables and ggplot2 figures displayed in the knitted report.
# ==============================================================================

prior_spec <- cfg$prior_specification
reduced_metrics <- reduced_comparison$metrics

software <- data.frame(
  component = c("R", "readxl", "ggplot2", "knitr", "rmarkdown"),
  version = c(R.version.string, as.character(utils::packageVersion("readxl")),
              as.character(utils::packageVersion("ggplot2")),
              as.character(utils::packageVersion("knitr")),
              as.character(utils::packageVersion("rmarkdown")))
)

configuration <- data.frame(
  setting = c("Loss threshold (percentage points)", "Rolling Nifty dates",
              "Minimum training exceedances", "Forecast step (Nifty dates)",
              "Cauchy slope scale", "Lasso rate", "Lasso smoothing epsilon",
              "Ridge precision", "Zellner g", "Numerical xi range", "Seed"),
  value = c(cfg$threshold, cfg$rolling_observations,
            cfg$minimum_training_exceedances, cfg$forecast_step,
            prior_spec$beta_cauchy_scale, prior_spec$beta_lasso_rate,
            prior_spec$beta_lasso_epsilon, prior_spec$beta_ridge_precision,
            prior_spec$g_prior_g,
            paste0("[", prior_spec$xi_lower, ", ", prior_spec$xi_upper, "]"),
            cfg$seed)
)

# 1. Data coverage check
merged <- main_merged  # [FIX B-03] restore the Date-typed merged data
data_checks_summary <- data.frame(
  Nifty_dates = nrow(merged),
  first_date = min(merged$Date),
  last_date = max(merged$Date),
  exceedances = sum(merged$exceedance == 1L, na.rm = TRUE),
  missing_nifty_return = sum(!is.finite(merged$nifty_return)),
  complete_prior_covariates = sum(complete.cases(merged[, c(
    "nifty_vol_lag1", "sp500_vol_lag1", "gold_vol_lag1", "india_vix_lag1", "usdinr_vol_lag1"
  )]))
)
stopifnot(!anyDuplicated(merged$Date), all(merged$excess[merged$exceedance == 1L] > 0))
knitr::kable(configuration, caption = "Shared analysis configuration")
knitr::kable(main_metrics, digits = 4, caption = "Rolling conditional severity scores")  # [FIX B-03]

# 2. Rolling prediction plot
if (nrow(main_predictions)) {
  p_rolling <- ggplot2::ggplot(main_predictions, ggplot2::aes(Date, predicted_mean_excess, color = model)) +
    ggplot2::geom_line(alpha = .7) +
    ggplot2::geom_point(ggplot2::aes(y = observed_excess), color = "black", alpha = .45, size = 1) +
    ggplot2::labs(x = NULL, y = "Observed excess / predicted conditional mean excess (pp)",  # [FIX B-04]
                  color = "Model", title = "Conditional severity forecasts on exceedance dates") +
    ggplot2::theme_minimal()
  print(p_rolling)
}

# 3. Rolling fit diagnostics summary
if (nrow(diagnostics)) {
  diag_summary <- do.call(rbind, lapply(split(diagnostics, diagnostics$model), function(d) {
    data.frame(model = d$model[1], attempted_fits = nrow(d),
               successful_fits = sum(!is.na(d$convergence) & d$convergence == 0L),
               failed_fits = sum(is.na(d$convergence) | d$convergence != 0L),
               min_training_exceedances = min(d$n_train))
  }))
  knitr::kable(diag_summary)
}

# 4. Reduced-model summary tables and plot
range_text <- function(x) {
  if (length(x) == 1L) return(sprintf("%.4f", x))
  sprintf("%.4f to %.4f", min(x, na.rm = TRUE), max(x, na.rm = TRUE))
}
set_summary <- do.call(rbind, lapply(split(reduced_metrics, reduced_metrics$model_set), function(d) {
  fit_rows <- reduced_comparison$diagnostics[
    reduced_comparison$diagnostics$model_set == d$model_set[1], , drop = FALSE
  ]
  vif_values <- fit_rows$max_vif[is.finite(fit_rows$max_vif)]
  data.frame(model_set = d$model_set[1], predictors = d$predictor_count[1],
             RMSE = range_text(d$RMSE), MAE = range_text(d$MAE),
             conditional_log_score = range_text(d$mean_conditional_log_score),
             max_rolling_VIF = if (length(vif_values)) round(max(vif_values), 3) else NA_real_,
             stringsAsFactors = FALSE)
}))
set_order <- c("constant_gpd", "full_5", "A", "B", "C", "D")
set_summary <- set_summary[match(set_order, set_summary$model_set), ]
knitr::kable(set_summary, caption = "Metric ranges across prior families by predictor set")
knitr::kable(reduced_metrics, digits = 4, caption = "Detailed results for all predictor/prior combinations")
knitr::kable(reduced_comparison$validation, caption = "Reduced-model comparison validation checks")

comparison_long <- do.call(rbind, lapply(c("RMSE", "MAE", "mean_conditional_log_score"), function(metric) {
  data.frame(model_set = reduced_metrics$model_set, prior = reduced_metrics$prior,
             metric = metric, value = reduced_metrics[[metric]], stringsAsFactors = FALSE)
}))
comparison_long$model_set <- factor(comparison_long$model_set, levels = set_order)
comparison_long$metric <- factor(comparison_long$metric,
  levels = c("RMSE", "MAE", "mean_conditional_log_score"),
  labels = c("RMSE (lower is better)", "MAE (lower is better)",
             "Conditional log score (higher is better)"))
p_reduced <- ggplot2::ggplot(comparison_long, ggplot2::aes(model_set, value, color = prior, group = prior)) +
  ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.55), size = 2.2) +
  ggplot2::facet_wrap(~metric, scales = "free_y", ncol = 1) +
  ggplot2::labs(x = "Predictor set", y = NULL, color = "Slope prior",
                title = "Reduced-model comparison on common rolling forecast dates") +
  ggplot2::theme_minimal()
print(p_reduced)

comparison_diag <- reduced_comparison$diagnostics
diag_summary_reduced <- do.call(rbind, lapply(split(comparison_diag,
                                           interaction(comparison_diag$model_set,
                                                       comparison_diag$prior, drop = TRUE)), function(d) {
  vif_values <- d$max_vif[is.finite(d$max_vif)]
  data.frame(model_set = d$model_set[1], prior = d$prior[1], attempted_fits = nrow(d),
             converged = sum(!is.na(d$convergence) & d$convergence == 0L),
             failed = sum(is.na(d$convergence) | d$convergence != 0L),
             median_max_vif = if (length(vif_values)) stats::median(vif_values) else NA_real_,
             max_vif = if (length(vif_values)) max(vif_values) else NA_real_)
}))
knitr::kable(diag_summary_reduced, digits = 3, caption = "Fit success and predictor collinearity by set and prior")
if (nrow(diagnostics)) knitr::kable(utils::head(diagnostics, 20), caption = "First 20 rolling fit diagnostics")
knitr::kable(validation, caption = "Runtime validation checks")
knitr::kable(software, caption = "Software versions used to knit this report")

# ==============================================================================
# STEP 10 (NEW): REPORT TABLES & FIGURES THAT WERE MISSING FROM THE PIPELINE
# Reproduces Tables 2-6, 8 and the conditional-probability values, and
# Figures 2-15 of the report, using only base R graphics (no extra packages).
# Inputs : results_dir/merged_data.csv and the reduced_model_comparison_*.csv files
# Outputs: results_dir/report_*.csv and report_dir/fig_*.png
# Figure 1 needs the raw OHLC price columns from the workbook and is not
# rebuilt here (see the audit log, bug B-15).
# ==============================================================================

report_dir <- file.path(output_base_dir, "report_figures")
dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)
rcsv <- function(f) utils::read.csv(file.path(results_dir, f), stringsAsFactors = FALSE)
wcsv <- function(x, f) utils::write.csv(x, file.path(results_dir, f), row.names = FALSE)
png_open <- function(f, w = 10, h = 5.5) grDevices::png(file.path(report_dir, f), width = w, height = h,
                                                         units = "in", res = 300, bg = "white")
cfg10 <- project_config(); ps10 <- cfg10$prior_specification; u10 <- cfg10$threshold
sets10 <- reduced_model_sets(); full10 <- sets10$full_5
md <- rcsv("merged_data.csv"); md$Date <- as.Date(md$Date)
ex10 <- md[md$exceedance == 1 & is.finite(md$excess) & md$excess > 0, ]
cc10 <- ex10[stats::complete.cases(ex10[, full10]), ]
rp <- rcsv("reduced_model_comparison_predictions.csv"); rp$Date <- as.Date(rp$Date)
rm10 <- rcsv("reduced_model_comparison_metrics.csv")
rw <- rcsv("reduced_model_comparison_windows.csv"); rw$forecast_date <- as.Date(rw$forecast_date)
pick <- function(set, prior) { d <- rp[rp$model_set == set & rp$prior == prior, ]; d[order(d$Date), ] }
col3 <- c(const = "#7F8C8D", C = "#2563EB", D = "#D97706", F5 = "#2E8B57")

# Unpenalised pooled GPD maximum-likelihood fit (Table 3/4). Note: the rolling
# "constant GPD" benchmark is a MAP fit with the Cauchy shape penalty; this is not.
gpd_mle <- function(z) {
  nll <- function(p) { v <- -sum(gpd_log_density(z, rep(p[1], length(z)), p[2])); if (is.finite(v)) v else 1e100 }
  f <- stats::nlminb(c(log(mean(z)), 0.1), nll, lower = c(-10, -5), upper = c(10, 1 - 1e-6))
  list(sigma = exp(f$par[1]), xi = f$par[2], loglik = -f$objective, convergence = f$convergence)
}
qgpd10 <- function(p, s, xi) if (abs(xi) < 1e-7) -s * log(1 - p) else s / xi * ((1 - p)^(-xi) - 1)
sgpd10 <- function(z, s, xi) if (abs(xi) < 1e-7) exp(-z / s) else pmax(1 + xi * z / s, 0)^(-1 / xi)

# ---- Table 2: descriptive statistics of excess losses
z10 <- ex10$excess
tab2 <- data.frame(statistic = c("Number of exceedances", "Mean excess", "Median excess", "SD", "Minimum",
                                 "Maximum", "Total losses > 5% (count)", "Total losses > 5% (share)"),
                   value = c(length(z10), mean(z10), median(z10), stats::sd(z10), min(z10), max(z10),
                             sum(z10 > 5 - u10), mean(z10 > 5 - u10)))
wcsv(tab2, "report_table2_descriptive.csv")

# ---- Table 3: pooled GPD (MLE)
pool <- gpd_mle(z10)
wcsv(data.frame(n = length(z10), xi = pool$xi, sigma = pool$sigma,
                implied_mean_excess = pool$sigma / (1 - pool$xi), loglik = pool$loglik),
     "report_table3_pooled_gpd.csv")

# ---- Table 4: threshold sensitivity
r10 <- md$nifty_return[is.finite(md$nifty_return)]
tab4 <- do.call(rbind, lapply(c(1.5, 1.75, 2, 2.25, 2.5, 2.75, 3), function(u) {
  zz <- -r10[r10 < -u] - u; f <- gpd_mle(zz)
  data.frame(threshold = u, exceedances = length(zz), mean_excess = mean(zz), xi_mle = f$xi, sigma_mle = f$sigma)
}))
wcsv(tab4, "report_table4_threshold_sensitivity.csv")

# ---- Table 8 (AIC/BIC) and Table 6 (Model C coefficients): full-sample MAP fits on complete cases
k_const <- fit_constant_gpd(cc10$excess, ps10)
ll_const <- sum(gpd_log_density(cc10$excess, rep(k_const$log_sigma, nrow(cc10)), k_const$xi))
ic <- list(data.frame(model_set = "constant_gpd", prior = "baseline", k = 2, loglik = ll_const))
coef_rows <- list()
for (s in names(sets10)) for (pr in c("cauchy", "lasso", "ridge", "g-prior")) {
  ft <- fit_gpd_map(cc10$excess, cc10[, sets10[[s]], drop = FALSE], pr, ps10)
  ll <- sum(gpd_log_density(cc10$excess, log(predict(ft, cc10[, sets10[[s]], drop = FALSE])), ft$xi))
  ic[[length(ic) + 1L]] <- data.frame(model_set = s, prior = pr, k = length(ft$beta) + 1L, loglik = ll)
  coef_rows[[length(coef_rows) + 1L]] <- data.frame(model_set = s, prior = pr,
    term = c("(Intercept)", sets10[[s]], "xi"), estimate = c(ft$beta, ft$xi))
}
ic <- do.call(rbind, ic); n10 <- nrow(cc10)
ic$AIC <- 2 * ic$k - 2 * ic$loglik; ic$BIC <- ic$k * log(n10) - 2 * ic$loglik
ic <- merge(ic, rm10[, c("model_set", "prior", "RMSE", "MAE", "mean_conditional_log_score")],
            by = c("model_set", "prior"), sort = FALSE)
wcsv(ic, "report_table8_aic_bic_rolling.csv")
coefs <- do.call(rbind, coef_rows); coefs$exp_estimate <- ifelse(coefs$term == "xi", NA, exp(coefs$estimate))
wcsv(coefs, "report_table6_full_sample_coefficients.csv")

# ---- Preferred model (chosen in Section 7.9 of the report): Model D, Cauchy prior.
# Conditional probabilities P(L>5 | L>2, x): India VIX at its 10th/50th/90th sample
# percentiles; USD/INR volatility held at its sample median. Change pref_* to re-run for another model.
pref_set <- "D"; pref_prior <- "cauchy"
pref_vars <- sets10[[pref_set]]; pref_label <- if (pref_vars[1] == "india_vix_lag1") "India VIX" else "Nifty volatility"
fitC <- fit_gpd_map(cc10$excess, cc10[, pref_vars], pref_prior, ps10)
std <- function(v, j) (v - fitC$means[j]) / fitC$sds[j]
zn <- std(cc10[[pref_vars[1]]], 1); zu <- std(cc10$usdinr_vol_lag1, 2)
probC <- function(a, b) { s <- exp(fitC$beta[1] + fitC$beta[2] * a + fitC$beta[3] * b); sgpd10(5 - u10, s, fitC$xi) }
qn <- stats::quantile(zn, c(.1, .5, .9))
wcsv(data.frame(model = paste(pref_set, pref_prior), percentile = c(10, 50, 90),
                first_covariate_raw = stats::quantile(cc10[[pref_vars[1]]], c(.1, .5, .9)),
                first_covariate_std = qn, usdinr_std_held_at = median(zu),
                P_L_gt5_given_L_gt2 = sapply(qn, probC, b = median(zu))),
     "report_conditional_probabilities_preferred.csv")

# ---- Figure 2: returns with exceedances
png_open("fig02_returns_exceedances.png", 11, 4.8)
par(mar = c(3, 4.5, 2.5, 1))
plot(md$Date, md$nifty_return, type = "l", col = "#7A8C9C", lwd = .5, xlab = "", ylab = "Daily log return (%)",
     main = sprintf("Nifty 50 daily log returns; %d exceedances below -%g%% marked", nrow(ex10), u10))
abline(h = -u10, col = "#C0392B", lty = 2, lwd = 1.5)
points(ex10$Date, ex10$nifty_return, pch = 19, cex = .45, col = "#C0392B")
dev.off()

# ---- Figure 3: histogram + pooled density
png_open("fig03_excess_histogram.png", 9, 5)
par(mar = c(4.5, 4.5, 2.5, 1))
hist(z10, breaks = seq(0, ceiling(max(z10)), by = 0.5), freq = FALSE, col = "#9EC1DE", border = "white",
     main = "Excess losses above the 2% threshold", xlab = "Excess loss z (percentage points)", ylab = "Density")
zg <- seq(0, max(z10), length.out = 400)
lines(zg, exp(gpd_log_density(pmax(zg, 1e-9), rep(log(pool$sigma), 400), pool$xi)), col = "#C0392B", lwd = 2)
legend("topright", bty = "n", lwd = c(10, 2), col = c("#9EC1DE", "#C0392B"),
       legend = c("Observed excesses (n = 179)", sprintf("Pooled GPD MLE (xi = %.3f, sigma = %.3f)", pool$xi, pool$sigma)))
dev.off()

# ---- Figure 4: survival (log scale)
png_open("fig04_survival.png", 9, 5)
par(mar = c(4.5, 4.5, 2.5, 1))
zs <- sort(z10); emp <- 1 - (seq_along(zs) - 0.5) / length(zs)
plot(zs, emp, log = "y", pch = 19, cex = .5, col = "#1F4E79", xlab = "Excess loss z (percentage points)",
     ylab = "P(Z > z)  [log scale]", main = "Empirical vs pooled GPD survival function", ylim = c(1e-3, 1))
lines(zg, sgpd10(zg, pool$sigma, pool$xi), col = "#C0392B", lwd = 2)
legend("topright", bty = "n", pch = c(19, NA), lty = c(NA, 1), lwd = c(NA, 2), col = c("#1F4E79", "#C0392B"),
       legend = c("Empirical (plotting position (i - 0.5)/n)", "Pooled GPD MLE"))
dev.off()

# ---- Figure 5: QQ plot
png_open("fig05_qq.png", 6.5, 6.2)
par(mar = c(4.5, 4.5, 2.5, 1))
th <- sapply((seq_along(zs) - 0.5) / length(zs), qgpd10, s = pool$sigma, xi = pool$xi)
plot(th, zs, pch = 19, cex = .6, col = "#1F4E79", xlab = "Theoretical pooled GPD quantile",
     ylab = "Empirical excess-loss quantile", main = "GPD Q-Q plot (pooled MLE)")
abline(0, 1, lty = 2, col = "#C0392B", lwd = 1.5)
dev.off()

# ---- Figure 6: threshold sensitivity
png_open("fig06_threshold_sensitivity.png", 12, 4)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 2.5, 1))
plot(tab4$threshold, tab4$exceedances, type = "b", pch = 19, xlab = "Threshold (%)", ylab = "Exceedances", main = "Exceedance count")
abline(v = u10, lty = 2, col = "grey40")
plot(tab4$threshold, tab4$mean_excess, type = "b", pch = 19, xlab = "Threshold (%)", ylab = "Mean excess (pp)", main = "Mean excess")
abline(v = u10, lty = 2, col = "grey40")
plot(tab4$threshold, tab4$xi_mle, type = "b", pch = 19, ylim = c(0, 0.35), xlab = "Threshold (%)",
     ylab = expression(hat(xi) ~ "(MLE)"), main = "Pooled shape estimate")
abline(v = u10, lty = 2, col = "grey40")
dev.off()

# ---- Figure 7: loss vs prior-date volatility (annualised by sqrt(250) for display)
png_open("fig07_loss_vs_volatility.png", 12, 4.3)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 2.5, 1))
loss <- -ex10$nifty_return
for (v in list(c("nifty_vol_lag1", "(a) Nifty 50"), c("sp500_vol_lag1", "(b) S&P 500"), c("gold_vol_lag1", "(c) Gold"))) {
  plot(ex10[[v[1]]] * sqrt(250), loss, pch = 19, cex = .6, col = "#2F6690",
       xlab = "Prior-date empirical volatility (annualised, %)", ylab = "Loss (-log return, %)", main = v[2])
  abline(h = c(2, 5), lty = 2, col = c("#D97706", "#C0392B"))
}
dev.off()

# ---- Figure 8: correlation heatmap on the 178 complete exceedances (the model's data)
labs10 <- c("Nifty vol", "S&P 500 vol", "Gold vol", "India VIX", "USD/INR vol")
cm <- stats::cor(cc10[, full10])
png_open("fig08_correlation_heatmap.png", 7, 6.3)
par(mar = c(6, 7, 3, 1))
pal <- grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(201)
image(1:5, 1:5, cm[, 5:1], col = pal, zlim = c(-1, 1), axes = FALSE, xlab = "", ylab = "",
      main = sprintf("Pearson correlations, %d complete exceedances", nrow(cc10)))
axis(1, 1:5, labs10, las = 2); axis(2, 1:5, rev(labs10), las = 1)
for (i in 1:5) for (j in 1:5) text(i, 6 - j, sprintf("%.2f", cm[i, j]), cex = .95)
dev.off()

# ---- Figure 9: training exceedances per window
png_open("fig09_training_exceedances.png", 10, 4.5)
par(mar = c(3, 4.5, 2.5, 1))
plot(rw$forecast_date, rw$n_training_exceedances, type = "b", pch = 19, col = "#1F4E79", xlab = "",
     ylab = "Training exceedances in window", main = "Training exceedances available at each forecast date")
abline(h = cfg10$minimum_training_exceedances, lty = 2, col = "#C0392B")
dev.off()

# ---- Figures 10-11: rolling shape and scale
three <- list(D = pick("D", "cauchy"), C = pick("C", "g-prior"), F5 = pick("full_5", "g-prior"))
lab3 <- c(D = "Model D (Cauchy)", C = "Model C (g-prior)", F5 = "Full-5 (g-prior)")
for (what in c("shape", "predicted_scale")) {
  png_open(if (what == "shape") "fig10_rolling_shape.png" else "fig11_rolling_scale.png", 10, 4.5)
  par(mar = c(3, 4.5, 2.5, 1))
  yl <- range(unlist(lapply(three, `[[`, what)))
  plot(three$D$Date, three$D[[what]], type = "n", ylim = yl, xlab = "",
       ylab = if (what == "shape") expression(hat(xi)) else expression(hat(sigma) ~ "(pp)"),
       main = if (what == "shape") "Rolling GPD shape estimates" else "Rolling predicted GPD scale")
  if (what == "shape") abline(h = 0, lty = 3)
  for (k in names(three)) points(three[[k]]$Date, three[[k]][[what]], pch = 19, cex = .8, col = col3[k])
  legend("topleft", bty = "n", pch = 19, col = col3[names(three)], legend = lab3)
  dev.off()
}

# ---- Figure 12: preferred model, full-sample fit
mC <- predict(fitC, cc10[, pref_vars]) / (1 - fitC$xi)
png_open("fig12_preferred_model_fit.png", 11, 4.8)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1))
plot(mC, cc10$excess, pch = 19, cex = .6, col = "#D97706", xlab = "Fitted conditional mean excess (pp)",
     ylab = "Observed excess (pp)", main = "(a) Observed vs fitted mean excess")
abline(0, 1, lty = 2)
grid_n <- seq(min(zn), max(zn), length.out = 200)
plot(grid_n, probC(grid_n, median(zu)), type = "l", lwd = 2, col = "#D97706", ylim = c(0, max(probC(grid_n, median(zu)))),
     xlab = paste("Standardised prior-date", pref_label), ylab = "P(L > 5% | L > 2%, x)",
     main = "(b) P(L > 5% | L > 2%, x)")
abline(v = qn, lty = 3, col = "grey40")
text(qn, 0, c("P10", "P50", "P90"), pos = 4, cex = .8)
dev.off()

# ---- Figure 13: observed vs predicted on the 28 dates
png_open("fig13_obs_vs_pred.png", 11, 4.8)
par(mar = c(3, 4.5, 2.5, 1))
cst <- pick("constant_gpd", "baseline")
plot(cst$Date, cst$observed_excess, pch = 19, col = "black", ylim = c(0, max(cst$observed_excess) * 1.05), xlab = "",
     ylab = "Excess loss (pp)", main = "Observed excess vs predicted conditional mean (28 forecast dates)")
points(cst$Date, cst$predicted_mean_excess, pch = 15, col = col3["const"])
points(three$C$Date, three$C$predicted_mean_excess, pch = 17, col = col3["C"])
points(three$D$Date, three$D$predicted_mean_excess, pch = 18, cex = 1.2, col = col3["D"])
legend("topleft", bty = "n", pch = c(19, 15, 17, 18), col = c("black", col3[c("const", "C", "D")]),
       legend = c("Observed", "Constant GPD", "Model C (g-prior)", "Model D (Cauchy)"))
dev.off()

# ---- Figure 14: metrics by set and prior (points, not lines across categories)
png_open("fig14_metrics_by_set.png", 12, 4.3)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 2.5, 1))
sets_order <- c("full_5", "A", "B", "C", "D"); pr4 <- c("cauchy", "lasso", "ridge", "g-prior")
pcol <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"); off <- c(-.15, -.05, .05, .15)
for (m in c("RMSE", "MAE", "mean_conditional_log_score")) {
  vals <- rm10[rm10$model_set %in% sets_order, ]
  base <- rm10[rm10$model_set == "constant_gpd", m]
  plot(NA, xlim = c(.6, 5.4), ylim = range(c(vals[[m]], base)), xaxt = "n", xlab = "Predictor set", ylab = m,
       main = if (m == "mean_conditional_log_score") "Mean log score (higher better)" else paste(m, "(lower better)"))
  axis(1, 1:5, c("Full-5", "A", "B", "C", "D")); abline(h = base, lty = 3)
  for (j in 1:4) { d <- vals[vals$prior == pr4[j], ]; points(match(d$model_set, sets_order) + off[j], d[[m]], pch = 19, col = pcol[j]) }
  if (m == "RMSE") legend("topright", bty = "n", pch = c(rep(19, 4), NA), lty = c(rep(NA, 4), 3),
                          col = c(pcol, "black"), legend = c(pr4, "Constant GPD"), cex = .85)
}
dev.off()

# ---- Figure 15: PIT diagnostics with KS p-values
pit_of <- function(d) 1 - mapply(sgpd10, d$observed_excess, d$predicted_scale, d$shape)
pits <- list(const = pit_of(cst), C = pit_of(three$C), D = pit_of(three$D))
ksp <- sapply(pits, function(p) suppressWarnings(stats::ks.test(p, "punif")$p.value))
wcsv(data.frame(model = names(pits), ks_p_value = ksp, median_pit = sapply(pits, median)), "report_pit_summary.csv")
png_open("fig15_pit.png", 11, 4.8)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1))
plot(c(0, 1), c(0, 1), type = "l", lty = 2, col = "grey40", xlab = "PIT value", ylab = "Empirical CDF", main = "PIT empirical CDF")
for (k in names(pits)) lines(stats::ecdf(pits[[k]]), col = col3[k], do.points = FALSE, verticals = TRUE, lwd = 1.6)
legend("topleft", bty = "n", lwd = 1.6, col = col3[names(pits)],
       legend = sprintf("%s (KS p = %.2f)", c("Constant GPD", "Model C (g-prior)", "Model D (Cauchy)"), ksp))
plot(c(0, 1), c(0, 1), type = "l", lty = 2, col = "grey40", xlab = "Uniform quantile", ylab = "Ordered PIT", main = "PIT Q-Q")
for (k in names(pits)) points((seq_along(pits[[k]]) - .5) / length(pits[[k]]), sort(pits[[k]]), pch = 19, cex = .6, col = col3[k])
dev.off()

message("Step 10 report tables written to ", results_dir, "; figures written to ", report_dir)

message("==============================================================================")
message("ALL STEPS COMPLETED SUCCESSFULLY!")
message("1. Core & Collinearity Results saved in: ", results_dir)
message("2. Presentation Tables & PNG Plots saved in: ", presentation_dir)
message("3. Report figures (Step 10) saved in: ", report_dir)
message("==============================================================================")
