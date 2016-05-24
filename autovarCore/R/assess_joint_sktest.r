#' Tests the skewness and kurtosis of a VAR model
#'
#' This function tests the joint skewness and kurtosis for the residuals of the endogenous variables in the specified VAR model. This function uses an implementation equivalent to STATA's \code{sktest}. Of the p-levels resulting from assessing the significance of the joint sktest for the residuals of that variable, the minimum is returned.
#' @param varest A \code{varest} model.
#' @return This function returns a p-level.
#' @examples
#' data_matrix <- matrix(nrow = 40, ncol = 3)
#' data_matrix[, ] <- runif(ncol(data_matrix) * nrow(data_matrix), 1, nrow(data_matrix))
#' colnames(data_matrix) <- c('rumination', 'happiness', 'activity')
#' varest <- autovarCore:::run_var(data_matrix, NULL, 1)
#' autovarCore:::assess_joint_sktest(varest)
assess_joint_sktest <- function(varest) {
  resids <- unname(resid(varest))
  nr_cols <- ncol(resids)
  nr_rows <- nrow(resids)
  if (is.null(nr_cols) || nr_cols < 1 || is.null(nr_rows) || nr_rows < 1)
    stop("No residuals found")
  minimum_p_level_sktest <- Inf
  for (column_index in 1:nr_cols) {
    column_resids <- resids[, column_index]
    coef_of_skewness <- coefficient_of_skewness(column_resids)
    coef_of_kurtosis <- coefficient_of_kurtosis(column_resids)
    z_skew <- z_skewness(coef_of_skewness, nr_rows)
    z_kurt <- z_kurtosis(coef_of_kurtosis, nr_rows)
    p_level_sktest <- sktest_joint_p(z_skew, z_kurt, nr_rows)
    if (p_level_sktest < minimum_p_level_sktest)
      minimum_p_level_sktest <- p_level_sktest
  }
  minimum_p_level_sktest
}

sktest_joint_p <- function(Z1, Z2, n) {
  K2 <- Z1 * Z1 + Z2 * Z2
  ZC2 <- -qnorm(exp(-0.5 * K2))
  logn <- log(n)
  cut <- 0.55 * (n^0.2) - 0.21
  a1 <- (-5 + 3.46 * logn) * exp(-1.37 * logn)
  b1 <- 1 + (0.854 - 0.148 * logn) * exp(-0.55 * logn)
  b2mb1 <- 2.13/(1 - 2.37 * logn)
  a2 <- a1 - b2mb1 * cut
  b2 <- b2mb1 + b1
  Z <- NULL
  if (ZC2 < -1) {
    Z <- ZC2
  } else if (ZC2 < cut) {
    Z <- a1 + b1 * ZC2
  } else {
    Z <- a2 + b2 * ZC2
  }
  P <- 1 - pnorm(Z)
  P
}
