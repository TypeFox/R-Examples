#' Wrapper function to estimate the error rate of a classifier
#'
#' We provide a wrapper function to estimate the error rate of a classifier
#' using any of the following estimators:#' 
#' \describe{
#'   \item{\code{\link{errorest_cv}}:}{Cross-validation Error Rate}
#'   \item{\code{\link{errorest_boot}}:}{Bootstrap Error Rate}
#'   \item{\code{\link{errorest_632}}:}{.632 Estimator from Efron (1983)}
#'   \item{\code{\link{errorest_632plus}}:}{.632+ Estimator from Efron and Tibshirani (1997)}
#'   \item{\code{\link{errorest_bcv}}:}{Bootstrap Cross-validation from Fu et al. (2005)}
#'   \item{\code{\link{errorest_loo_boot}}:}{Leave-One-Out Bootstrap Error Rate}
#'   \item{\code{\link{errorest_apparent}}:}{Apparent Error Rate}
#' }
#'
#' This wrapper function provides a common means to estimate classification
#' error rates and is useful for simulation studies where multiple error-rate
#' estimators are being considered.
#' 
#' For details about an individual error-rate estimator, please see its
#' respective documentation.
#'
#' @param x a matrix of n observations and p features
#' @param y a vector of n class labels. (Must to be a 'factor'.)
#' @param estimator the estimator used to compute the error rate
#' @param train a function that builds the classifier. (See details.)
#' @param classify a function that classifies observations from the constructed
#' classifier from \code{train}. (See Details.)
#' @param ... additional arguments passed to the error-rate estimation code
#' @return an estimate of the classifier's error rate
#' 
#' @import MASS
#' @export
#' @examples
#' require('MASS')
#' iris_x <- data.matrix(iris[, -5])
#' iris_y <- iris[, 5]
#'
#' # Because the \code{classify} function returns multiples objects in a list,
#' # we provide a wrapper function that returns only the class labels.
#' lda_wrapper <- function(object, newdata) { predict(object, newdata)$class }
#'
#' # Cross-Validation (default)
#' errorest(x = iris_x, y = iris_y, train = MASS:::lda, classify = lda_wrapper)
#'
#' # .632
#' errorest(x = iris_x, y = iris_y, estimator = "632", train = MASS:::lda,
#'          classify = lda_wrapper)
#'
#' # Bootstrap Error Rate
#' # The argument 'num_bootstraps' is passed on to 'errorest_boot'
#' errorest(x = iris_x, y = iris_y, estimator = "boot", train = MASS:::lda,
#'          classify = lda_wrapper, num_bootstraps = 42)
#' 
errorest <- function(x, y, estimator = c("cv", "boot", "632", "632+", "bcv",
                               "loo-boot", "apparent"), train, classify, ...) {
  x <- as.matrix(x)
  y <- as.factor(y)

  estimator <- match.arg(estimator)

  switch(estimator,
    cv = errorest_cv(x = x, y = y, train = train, classify = classify, ...),
    boot = errorest_boot(x = x, y = y, train = train, classify = classify, ...),
    `632` = errorest_632(x = x, y = y, train = train, classify = classify, ...),
    `632+` = errorest_632plus(x = x, y = y, train = train, classify = classify, ...),
    bcv = errorest_bcv(x = x, y = y, train = train, classify = classify, ...),
    `loo-boot` = errorest_loo_boot(x = x, y = y, train = train, classify = classify, ...),
    apparent = errorest_apparent(x = x, y = y, train = train, classify = classify, ...)
  )
}
