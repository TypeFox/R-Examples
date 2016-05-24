#' Return the model fit for the given varest model
#'
#' This function returns the model fit for the given model as either an AIC or BIC score. We compensating for logtransformation so that the model scores of logtransformed and non-logtransformed models can be compared with each other directly. This compensation is implemented by subtracting the logtransformed data from the log-likelihood score and using the result as log-likelihood score for the AIC/BIC calculations.
#' @param varest A \code{varest} model.
#' @param criterion A character string being either \code{'AIC'} or \code{'BIC'}.
#' @param logtransformed A boolean, either \code{TRUE} or \code{FALSE}, indicating whether the input data for the model has been logtransformed.
#' @return This returns a floating point that is either the AIC or BIC criterion for the model. A lower number corresponds to a better model fit.
#' @examples
#' data_matrix <- matrix(nrow = 40, ncol = 3)
#' data_matrix[, ] <- runif(ncol(data_matrix) * nrow(data_matrix), 1, nrow(data_matrix))
#' colnames(data_matrix) <- c('rumination', 'happiness', 'activity')
#' varest <- autovarCore:::run_var(data_matrix, NULL, 1)
#' autovarCore:::model_score(varest, 'AIC', FALSE)
#' @export
model_score <- function(varest, criterion, logtransformed) {
  if (!(criterion %in% c('AIC', 'BIC')))
    stop(paste("Unknown criterion:", criterion))
  nr_observations <- summary(varest)$obs
  nr_estimated_params <- nr_estimated_parameters(varest)
  loglikelihood <- determine_loglikelihood(varest, logtransformed)
  switch(criterion,
         'AIC' = -2 * loglikelihood + 2 * nr_estimated_params,
         'BIC' = -2 * loglikelihood + log(nr_observations) * nr_estimated_params)
}

determine_loglikelihood <- function(varest, logtransformed) {
  if (logtransformed)
    loglikelihood_for_logtransformed(varest)
  else
    summary(varest)$logLik
}

loglikelihood_for_logtransformed <- function(varest) {
  # http://webspace.qmul.ac.uk/aferreira/lect2-var2_handout.pdf
  # http://www.unc.edu/courses/2010fall/ecol/563/001/docs/lectures/lecture15.htm#transformation
  nr_observations <- varest$obs
  var_dimension <- varest$K
  resids <- resid(varest)
  sigma <- crossprod(resids)/nr_observations
  r <- -(nr_observations * var_dimension/2) * log(2 * pi) - (nr_observations/2) * log(det(sigma)) -
    (1/2) * sum(diag(resids %*% solve(sigma) %*% t(resids)))
  r <- r - sum(varest$y)
  r
}

nr_estimated_parameters <- function(varest) {
  result <- 0
  for (equation in varest$varresult)
    result <- result + length(equation$coefficients)
  result
}
