#' Tests the homeskedasticity assumption for a VAR model using a portmanteau test on the squared residuals
#'
#' This function tests the homeskedasticity assumption for the residuals of the endogenous variables in the specified VAR model. This function implements the portmanteau squared test known as the Ljung-Box test, and results are comparable with STATA's \code{wntestq}. Of the p-levels resulting from assessing the homeskedasticity assumption for the squared residuals of that variable, the minimum is returned.
#' @param varest A \code{varest} model.
#' @return This function returns a p-level.
#' @examples
#' data_matrix <- matrix(nrow = 40, ncol = 3)
#' data_matrix[, ] <- runif(ncol(data_matrix) * nrow(data_matrix), 1, nrow(data_matrix))
#' colnames(data_matrix) <- c('rumination', 'happiness', 'activity')
#' varest <- autovarCore:::run_var(data_matrix, NULL, 1)
#' autovarCore:::assess_portmanteau_squared(varest)
assess_portmanteau_squared <- function(varest) {
  data <- unname(resid(varest))
  portmanteau_test_data(data^2)
}
