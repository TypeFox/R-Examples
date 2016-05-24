#' Execute a series of model validity assumptions
#'
#' This function returns the given suite of tests on for the given VAR model. For each test, the result is the minimum p-level of all the assumptions and p-levels checked within the test. In other words, the result of a test is the p-level that should be used as a threshold below which outcomes are considered statistically significant (e.g., a result of 0.06 is better than a result of 0.03). The \code{run_tests} function returns a vector of results, one for each test, in the order corresponding to the \code{test_names} argument.
#' @param varest A \code{varest} model.
#' @param test_names A vector of names of tests given as character strings. Supported tests are specified in the \code{autovarCore:::supported_test_names()} function.
#' @return This function returns a vector of p-levels.
#' @examples
#' data_matrix <- matrix(nrow = 40, ncol = 3)
#' data_matrix[, ] <- runif(ncol(data_matrix) * nrow(data_matrix), 1, nrow(data_matrix))
#' colnames(data_matrix) <- c('rumination', 'happiness', 'activity')
#' varest <- autovarCore:::run_var(data_matrix, NULL, 1)
#' autovarCore:::run_tests(varest, 'portmanteau')
#' @export
run_tests <- function(varest, test_names) {
  results <- NULL
  for (test_name in test_names)
    results <- c(results, run_test(test_name)(varest))
  results
}
