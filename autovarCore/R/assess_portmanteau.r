#' Tests the white noise assumption for a VAR model using a portmanteau test on the residuals
#'
#' This function tests the white noise assumption for the residuals of the endogenous variables in the specified VAR model. This function implements the portmanteau test known as the Ljung-Box test, and results are comparable with STATA's \code{wntestq}. Of the p-levels resulting from assessing the white noise assumption for the residuals of that variable, the minimum is returned.
#' @param varest A \code{varest} model.
#' @return This function returns a p-level.
#' @examples
#' data_matrix <- matrix(nrow = 40, ncol = 3)
#' data_matrix[, ] <- runif(ncol(data_matrix) * nrow(data_matrix), 1, nrow(data_matrix))
#' colnames(data_matrix) <- c('rumination', 'happiness', 'activity')
#' varest <- autovarCore:::run_var(data_matrix, NULL, 1)
#' autovarCore:::assess_portmanteau(varest)
assess_portmanteau <- function(varest) {
  data <- unname(resid(varest))
  portmanteau_test_data(data)
}

portmanteau_test_data <- function(data) {
  # This function is also used by assess_portmanteau_squared.
  nr_cols <- ncol(data)
  nr_rows <- nrow(data)
  if (is.null(nr_cols) || nr_cols < 1 || is.null(nr_rows) || nr_rows < 1)
    stop("No residuals found")
  port_lags <- determine_portmanteau_lags(data)
  if (port_lags < 1)
    stop("Not enough observations in the data")
  minimum_p_level_port <- Inf
  for (column_index in 1:nr_cols) {
    column_data <- data[, column_index]
    port_test_statistic <- portmanteau_test_statistic(column_data, nr_rows, port_lags)
    p_level_port <- chi_squared_prob(port_test_statistic, port_lags)
    if (p_level_port < minimum_p_level_port)
      minimum_p_level_port <- p_level_port
  }
  minimum_p_level_port
}

determine_portmanteau_lags <- function(data) {
  # This is the default value used in STATA.
  min(floor(nrow(data)/2) - 2, 40)
}

portmanteau_test_statistic <- function(data, n, h) {
  data <- data - mean(data)
  suma <- 0
  for (k in 1:h)
    suma <- suma + (sample_autocorrelation(data, k, n)^2)/(n - k)
  q <- n * (n + 2) * suma
  q
}

sample_autocorrelation <- function(data, k, n) {
  res <- 0
  for (t in (k + 1):n)
    res <- res + data[t] * data[t - k]
  # See the paper of Ljung-Box test for this definition of autocorrelation.
  denom <- 0
  for (t in 1:n)
    denom <- denom + data[t]^2
  res <- res/denom
  res
}

chi_squared_prob <- function(q, h) {
  pchisq(q, h, lower.tail = FALSE)
}
