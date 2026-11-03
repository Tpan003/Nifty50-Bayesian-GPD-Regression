# Data import, date checks, lagged predictor construction, and exceedances.

read_market_sheets <- function(path) {
  expected_sheets <- c(
    "NIFTY 50 DATA INR",
    "S&P500 DATA USD",
    "GOLD DATA USD",
    "INDIA VIX DATA INR",
    "USD-INR DATA"
  )

  available_sheets <- readxl::excel_sheets(path)
  missing_sheets <- setdiff(expected_sheets, available_sheets)
  if (length(missing_sheets) > 0L) {
    stop("Workbook is missing sheets: ", paste(missing_sheets, collapse = ", "))
  }

  read_one_sheet <- function(sheet_name) {
    # A period is used as a missing-value marker in some workbook cells.
    suppressWarnings(
      readxl::read_excel(path, sheet = sheet_name, na = c("", "."))
    )
  }

  sheets <- lapply(expected_sheets, read_one_sheet)
  names(sheets) <- c("nifty", "sp500", "gold", "vix", "usdinr")
  sheets
}

check_dates <- function(data, label = "Dataset") {
  if (!"Date" %in% names(data)) {
    stop(label, " does not contain a Date column.")
  }

  data$Date <- as.Date(data$Date)
  if (anyNA(data$Date)) {
    stop(label, " contains missing or invalid dates.")
  }
  if (anyDuplicated(data$Date)) {
    stop(label, " contains duplicate dates.")
  }

  data[order(data$Date), , drop = FALSE]
}

# For a forecast made on date t, return the latest observation dated strictly
# before t. This allows markets with different holiday calendars to be joined.
latest_prior <- function(target_dates, source_dates, values) {
  stopifnot(length(source_dates) == length(values))

  valid <- !is.na(source_dates) & is.finite(values)
  source_dates <- as.Date(source_dates[valid])
  values <- values[valid]

  order_index <- order(source_dates)
  source_dates <- source_dates[order_index]
  values <- values[order_index]

  target_dates <- as.Date(target_dates)
  index <- findInterval(
    as.numeric(target_dates - 1),
    as.numeric(source_dates)
  )

  result <- rep(NA_real_, length(target_dates))
  available <- index > 0L
  result[available] <- values[index[available]]
  result
}

prepare_and_merge <- function(path, threshold = 2) {
  if (length(threshold) != 1L || !is.finite(threshold) || threshold <= 0) {
    stop("threshold must be one positive finite number.")
  }

  sheets <- read_market_sheets(path)
  nifty_sheet <- sheets$nifty

  required_nifty_columns <- c(
    "DATE", "Log return (%)", "Daily empirical volatility (%)"
  )
  missing_columns <- setdiff(required_nifty_columns, names(nifty_sheet))
  if (length(missing_columns) > 0L) {
    stop("Nifty sheet is missing columns: ", paste(missing_columns, collapse = ", "))
  }

  nifty <- data.frame(
    Date = as.Date(nifty_sheet$DATE),
    nifty_return = as.numeric(nifty_sheet[["Log return (%)"]]),
    nifty_vol = as.numeric(nifty_sheet[["Daily empirical volatility (%)"]])
  )
  nifty <- check_dates(nifty, "Nifty sheet")

  # Nifty volatility is lagged as well. The function excludes same-day values.
  nifty$nifty_vol_lag1 <- latest_prior(
    nifty$Date, nifty$Date, nifty$nifty_vol
  )

  predictor_sources <- list(
    sp500 = list(
      sheet = sheets$sp500,
      value_column = "Daily empirical volatility (%)",
      output_column = "sp500_vol_lag1"
    ),
    gold = list(
      sheet = sheets$gold,
      value_column = "Daily empirical volatility (%)",
      output_column = "gold_vol_lag1"
    ),
    vix = list(
      sheet = sheets$vix,
      value_column = "India VIX level (points)",
      output_column = "india_vix_lag1"
    ),
    usdinr = list(
      sheet = sheets$usdinr,
      value_column = "Daily empirical volatility (%)",
      output_column = "usdinr_vol_lag1"
    )
  )

  for (source_name in names(predictor_sources)) {
    source_info <- predictor_sources[[source_name]]
    source_sheet <- source_info$sheet
    date_column <- if ("DATE" %in% names(source_sheet)) "DATE" else "Date"

    required_columns <- c(date_column, source_info$value_column)
    missing_columns <- setdiff(required_columns, names(source_sheet))
    if (length(missing_columns) > 0L) {
      stop(
        "The ", source_name, " sheet is missing columns: ",
        paste(missing_columns, collapse = ", ")
      )
    }

    source_data <- data.frame(
      Date = as.Date(source_sheet[[date_column]]),
      value = as.numeric(source_sheet[[source_info$value_column]])
    )
    source_data <- check_dates(source_data, source_name)

    nifty[[source_info$output_column]] <- latest_prior(
      nifty$Date, source_data$Date, source_data$value
    )
  }

  nifty$threshold <- threshold
  nifty$exceedance <- as.integer(
    is.finite(nifty$nifty_return) & nifty$nifty_return < -threshold
  )
  nifty$excess <- ifelse(
    nifty$exceedance == 1L,
    -nifty$nifty_return - threshold,
    NA_real_
  )

  # If the workbook already contains an excess-loss column, compare it with
  # the value calculated from returns. The calculated value remains primary.
  stored_column <- paste0(
    "Excess loss beyond ", threshold, "% threshold (pp)"
  )
  if (stored_column %in% names(nifty_sheet)) {
    stored <- data.frame(
      Date = as.Date(nifty_sheet$DATE),
      stored_excess = as.numeric(nifty_sheet[[stored_column]])
    )
    stored <- check_dates(stored, "Stored excess-loss column")
    stored_values <- stored$stored_excess[match(nifty$Date, stored$Date)]
    comparable <- is.finite(stored_values) & is.finite(nifty$excess)

    if (any(comparable) &&
        max(abs(stored_values[comparable] - nifty$excess[comparable])) > 1e-6) {
      warning(
        "The stored excess-loss column differs from the recomputed values. ",
        "The recomputed values will be used."
      )
    }
  }

  nifty
}

build_extremes <- function(merged_data) {
  extremes <- merged_data[
    merged_data$exceedance == 1L &
      is.finite(merged_data$excess) &
      merged_data$excess > 0,
    , drop = FALSE
  ]

  extremes <- extremes[order(extremes$Date), , drop = FALSE]
  extremes
}
