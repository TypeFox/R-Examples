#' Tests the kurtosis of a VAR model
#'
#' This function tests the kurtosis for the residuals of the endogenous variables in the specified VAR model. This function uses an implementation equivalent to STATA's \code{sktest}. Of the p-levels resulting from assessing the significance of the kurtosis for the residuals of that variable, the minimum is returned.
#' @param varest A \code{varest} model.
#' @return This function returns a p-level.
#' @examples
#' data_matrix <- matrix(nrow = 40, ncol = 3)
#' data_matrix[, ] <- runif(ncol(data_matrix) * nrow(data_matrix), 1, nrow(data_matrix))
#' colnames(data_matrix) <- c('rumination', 'happiness', 'activity')
#' varest <- autovarCore:::run_var(data_matrix, NULL, 1)
#' autovarCore:::assess_kurtosis(varest)
assess_kurtosis <- function(varest) {
  resids <- unname(resid(varest))
  nr_cols <- ncol(resids)
  nr_rows <- nrow(resids)
  if (is.null(nr_cols) || nr_cols < 1 || is.null(nr_rows) || nr_rows < 1)
    stop("No residuals found")
  minimum_p_level_kurt <- Inf
  for (column_index in 1:nr_cols) {
    column_resids <- resids[, column_index]
    coef_of_kurtosis <- coefficient_of_kurtosis(column_resids)
    z_kurt <- z_kurtosis(coef_of_kurtosis, nr_rows)
    p_level_kurt <- 2 - 2 * pnorm(abs(z_kurt))
    if (p_level_kurt < minimum_p_level_kurt)
      minimum_p_level_kurt <- p_level_kurt
  }
  minimum_p_level_kurt
}

z_kurtosis <- function(b2, n) {
  # This function is also used by assess_joint_sktest.
  Eb2 <- (3 * (n - 1)/(n + 1))
  varb2 <- (24 * n * (n - 2) * (n - 3))/((n + 1)^2 * (n + 3)*(n + 5))
  X <- (b2 - Eb2) / sqrt(varb2)
  sqrtB1b2 <- (6 * (n^2 - 5 * n + 2)/((n + 7) * (n + 9))) *
              sqrt((6 * (n + 3) * (n + 5))/(n * (n - 2) * (n - 3)))
  A <- 6 + (8/sqrtB1b2) * ((2/sqrtB1b2) + sqrt(1 + (4/(sqrtB1b2^2))))
  Z2 <- (1/sqrt(2/(9 * A))) * ((1 - (2/(9 * A))) - ((1 - (2/A))/(1 + X * sqrt(2/(A - 4))))^(1/3))
  Z2
}

coefficient_of_kurtosis <- function(x) {
  # This function is also used by assess_joint_sktest.
  m4 <- rth_moment_about_the_mean(x, 4)
  m2 <- rth_moment_about_the_mean(x, 2)
  m4 * m2^(-2)
}
