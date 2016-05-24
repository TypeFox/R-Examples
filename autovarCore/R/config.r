# Configuration and defaults

default_autovar_params <- function() {
  list(significance_levels = c(0.05, 0.01, 0.005),
       test_names = c('portmanteau', 'portmanteau_squared', 'skewness'),
       criterion = 'AIC',
       imputation_iterations = 30,
       measurements_per_day = 1)
}

supported_test_names <- function() {
  c('portmanteau',
    'portmanteau_squared',
    'skewness',
    'kurtosis',
    'joint_sktest')
}

run_test <- function(test_name) {
  test_function <- switch(test_name,
         'portmanteau' = assess_portmanteau,
         'portmanteau_squared' = assess_portmanteau_squared,
         'skewness' = assess_skewness,
         'kurtosis' = assess_kurtosis,
         'joint_sktest' = assess_joint_sktest)
  if (is.null(test_function))
    stop(paste("Unknown test:", test_name))
  test_function
}

supported_criteria <- function() {
  c('AIC',
    'BIC')
}

p_level_for_trend_significance <- function() {
  0.05
}

std_factor_for_normal_outliers <- function() {
  2.5
}

std_factor_for_squared_outliers <- function() {
  2.5
}
