# Forecasting Extreme Nifty 50 Loss Severity

A reproducible R project for modelling the severity of extreme Nifty 50 daily losses using a Peaks-Over-Threshold (POT) approach and the Generalised Pareto Distribution (GPD).

The code has been organised into separate files so that the data preparation, model estimation, validation, reduced-model comparison, diagnostics, and reporting steps are easier to read and maintain.

## What the project does

1. Imports Nifty 50, S&P 500, gold, India VIX, and USD/INR data from one Excel workbook.
2. Checks the date columns and aligns each predictor to the latest observation strictly before the Nifty forecast date.
3. Defines an extreme-loss event when the Nifty daily log return is below `-2` percentage points. The excess response is the loss beyond that threshold.
4. Fits an intercept-only GPD benchmark and covariate-dependent GPD scale models using Cauchy, lasso, ridge, and Zellner g-prior regularisation.
5. Evaluates the models using rolling training windows and held-out exceedance dates.
6. Compares the full five-predictor model with four prespecified reduced predictor sets.
7. Runs data, timing, likelihood, fit, and multicollinearity diagnostics.
8. Produces CSV tables, presentation charts, report figures, and an optional knitted HTML report.

**Important:** the model forecasts excess-loss severity conditional on a threshold exceedance. It does not forecast the probability of an exceedance occurring.

## Repository structure

```text
nifty50-extreme-losses/
├── R/
│   ├── data_functions.R
│   ├── model_functions.R
│   ├── project_config.R
│   ├── validation_functions.R
│   ├── reduced_model_comparison.R
│   ├── collinearity_diagnostics.R
│   ├── presentation_materials.R
│   └── report_materials.R
├── data/
│   └── README.md
├── results/                 # Generated output; not committed to Git
├── presentation_materials/  # Generated output; not committed to Git
├── report_figures/          # Generated output; not committed to Git
├── run_analysis.R
├── report.Rmd
├── .gitignore
└── README.md
```

## Requirements

Use R 4.1 or later. The project uses these packages:

- `readxl` for importing the Excel workbook.
- `ggplot2` for presentation charts.

The report can optionally be knitted to HTML using `rmarkdown` and `knitr`.

Install the required packages:

```r
install.packages(c("readxl", "ggplot2"))
```

For the optional HTML report:

```r
install.packages(c("rmarkdown", "knitr"))
```

## Input workbook

Place the workbook in the `data/` folder, or provide its full path when running the script. The expected worksheet names are:

| Worksheet | Required columns |
|---|---|
| `NIFTY 50 DATA INR` | `DATE`, `Log return (%)`, `Daily empirical volatility (%)` |
| `S&P500 DATA USD` | `DATE` or `Date`, `Daily empirical volatility (%)` |
| `GOLD DATA USD` | `DATE` or `Date`, `Daily empirical volatility (%)` |
| `INDIA VIX DATA INR` | `DATE` or `Date`, `India VIX level (points)` |
| `USD-INR DATA` | `DATE` or `Date`, `Daily empirical volatility (%)` |

The date columns should be recognised by Excel/R as dates. Empty cells and cells containing a period (`.`) are treated as missing values.

The workbook itself is not included in this repository. Add your authorised project dataset locally. Do not upload private, restricted, or unlicensed market data to a public repository.

## Run the analysis

Open the project folder as your working directory in RStudio, then run:

```r
source("run_analysis.R")
```

In an interactive session, the script asks you to select the input workbook and output folder. On systems without a folder-selection dialog, it uses the current working directory for outputs.

You can also run it from a terminal or the RStudio Terminal:

```bash
Rscript run_analysis.R "data/RT_ONLY_COVARIATES.xlsx" "."
```

The first argument is the Excel workbook. The second argument is the output folder. If the second argument is omitted, results are written to the current working directory.

Environment variables are also supported:

```bash
NIFTY_INPUT="data/RT_ONLY_COVARIATES.xlsx" NIFTY_OUTPUT="." Rscript run_analysis.R
```

## Analysis settings

The main settings are defined in `R/project_config.R`:

| Setting | Default |
|---|---:|
| Extreme-loss threshold | 2 percentage points |
| Rolling window | 1,250 prior Nifty trading dates |
| Minimum training exceedances | 40 |
| Forecast step | 1 Nifty trading date |
| Random seed | 20261018 |

The predictor columns are lagged using the latest available source observation strictly before the Nifty forecast date. This is important because market calendars do not always match. For example, a US market observation dated on the Nifty forecast date is not used for that forecast.

## Model specification

For a threshold exceedance, the excess loss is modelled using a GPD with scale depending on lagged predictors:

\[
\log(\sigma_t) = \beta_0 + \beta_1 x_{1,t-1} + \cdots + \beta_p x_{p,t-1}.
\]

The shape parameter is estimated subject to numerical bounds. Predictor columns are standardised using the training sample only. The training means and standard deviations are then reused for the held-out prediction.

The slope regularisation options are:

- Cauchy prior penalty
- Smoothed lasso penalty
- Ridge penalty
- Zellner g-prior penalty

The reduced predictor specifications are prespecified in `R/reduced_model_comparison.R`:

| Set | Predictors |
|---|---|
| Full 5 | Nifty volatility, S&P 500 volatility, gold volatility, India VIX, USD/INR volatility |
| A | Nifty volatility, S&P 500 volatility, gold volatility, USD/INR volatility |
| B | S&P 500 volatility, gold volatility, India VIX, USD/INR volatility |
| C | Nifty volatility, USD/INR volatility |
| D | India VIX, USD/INR volatility |

Each set is evaluated with all four prior specifications on the same rolling forecast dates. The intercept-only GPD is included as a benchmark.

## Output files

The script creates three output folders inside the chosen output directory.

### `results/`

Contains the main model outputs and diagnostics, including:

- `merged_data.csv`
- `extreme_severity_data.csv`
- `rolling_predictions.csv`
- `rolling_metrics.csv`
- `fit_diagnostics.csv`
- `validation_checks.csv`
- `reduced_model_comparison_predictions.csv`
- `reduced_model_comparison_metrics.csv`
- `reduced_model_comparison_fit_diagnostics.csv`
- `reduced_model_comparison_windows.csv`
- `reduced_model_comparison_validation.csv`
- `collinearity_*.csv` and `collinearity_report.md`
- `report_table*.csv` and `report_*.csv`
- `analysis_objects.rds` and `session_info.txt`

### `presentation_materials/`

Contains presentation-ready CSV tables, five PNG charts, and a short `README.md` describing the chart and metric definitions.

### `report_figures/`

Contains the additional base-R report figures (`fig02` through `fig15`). Figure 1 is not reconstructed by this pipeline because the supplied input specification does not require the raw OHLC price columns used for it.

## Create the HTML report

After the analysis has completed, knit the report from the repository root:

```r
rmarkdown::render(
  "report.Rmd",
  params = list(output_dir = ".")
)
```

Change `output_dir` if you saved the generated folders somewhere else. The report includes the rolling metrics, validation results, presentation charts, and additional report figures.

## Interpreting the results

- **RMSE and MAE** compare the predicted conditional mean excess loss with the observed excess loss on held-out exceedance dates.
- **Conditional log score** evaluates the fitted GPD density for the observed excess. A higher (less negative) score is better.
- **Out-of-sample R-squared** is calculated as `1 - SSE/SST` over the forecast dates. It can be negative and is not adjusted R-squared.
- **VIF and pairwise correlation** describe collinearity among predictors. They do not, on their own, prove overfitting.
- **AIC/BIC report table** uses full-sample MAP parameter estimates with the corresponding likelihood and parameter count, then places the rolling metrics beside them. Treat these as supplementary comparisons, not as a substitute for out-of-sample evaluation.
- **PIT and KS diagnostics** are based on a relatively small set of held-out exceedance forecasts. Interpret them cautiously.

The reduced-model comparison dates are also used to compare specifications. A model selected using those results should be treated as exploratory until it is tested on a separate, untouched period or within a nested rolling validation procedure.

## Reproducibility notes

- The project configuration and prior hyperparameters are stored with the saved analysis objects.
- The run saves `sessionInfo()` output so package and R versions can be checked later.
- The random seed is set from the shared project configuration.
- Dates are exported to CSV as ISO-formatted `YYYY-MM-DD` strings.

## Current scope note

The supplied pipeline defines report outputs for Tables 2, 3, 4, 6, and 8. A calculation or specification for Table 5 was not present in the source code, so this repository does not invent one. Add the intended Table 5 definition to `R/report_materials.R` when its contents are confirmed.

## GitHub upload

1. Create a new repository on GitHub.
2. Unzip the project folder locally.
3. Open the extracted folder in RStudio or File Explorer.
4. Upload the project files and folders. Do not upload the raw workbook or generated outputs unless you intend to publish them.
5. Commit with a message such as `Add reproducible Nifty 50 extreme-loss analysis pipeline`.

