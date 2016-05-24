#' Determines if a trend is required for the specified VAR model
#'
#' This function uses the Phillips-Perron Unit Root Test to determine whether a trend is required for a VAR model based on the given matrix of endogenous variables and the given lag. All variables are assessed individually. This function returns \code{TRUE} if any of the endogenous variables requires a trend.
#' @param endo_matrix The matrix of endogenous variables in the model.
#' @param lag An integer specifying the lag length of the model.
#' @return A boolean indicating whether a trend is required for the specified VAR model.
#' @examples
#' data_matrix <- matrix(nrow = 40, ncol = 3)
#' data_matrix[, ] <- runif(ncol(data_matrix) * nrow(data_matrix), 1, 10)
#' data_matrix[, 3] <- (1:40) + rnorm(40)
#' colnames(data_matrix) <- c('rumination', 'happiness', 'activity')
#' data_matrix
#' autovarCore:::needs_trend(data_matrix, 1)
#' @export
needs_trend <- function(endo_matrix, lag) {
  for (column_index in 1:ncol(endo_matrix))
    if (column_needs_trend(endo_matrix[, column_index], lag))
      return(TRUE)
  FALSE
}

column_needs_trend <- function(endo_column, lag) {
  test_type <- "Z-tau"
  deterministic_part_in_test_regression <- "trend"
  pp_summary <- urca::summary(ur.pp(x = endo_column,
                                    type = test_type,
                                    model = deterministic_part_in_test_regression,
                                    use.lag = lag))
  trend_p <- pp_summary@testreg$coefficients[rownames(pp_summary@testreg$coefficients) == 'trend', ][[4]]
  if (trend_p > p_level_for_trend_significance())
    return(FALSE)
  five_crit_val <- pp_summary@cval[colnames(pp_summary@cval) == '5pct']
  teststat <- pp_summary@teststat
  # If the p-level of trend is significant AND the Z(t) value is less than
  # or equal to the 5% critical value, then the variable requires a trend.
  teststat <= five_crit_val
}
