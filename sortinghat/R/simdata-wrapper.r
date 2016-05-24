#' Wrapper function to generate data from a variety of data-generating families
#' for classification studies.
#'
#' We provide a wrapper function to generate random variates from any of the
#' following data-generating families:
#' \describe{
#'   \item{\code{\link{simdata_normal}}:}{Multivariate normal}
#'   \item{\code{\link{simdata_t}}:}{Multivariate Student's t}
#'   \item{\code{\link{simdata_uniform}}:}{Multivariate uniform}
#'   \item{\code{\link{simdata_contaminated}}:}{Multivariate contaminated normal}
#'   \item{\code{\link{simdata_guo}}:}{Simulation configuration from Guo et al. (2007)}
#'   \item{\code{\link{simdata_friedman}}:}{Six simulation configurations from Friedman (1989)}
#' }
#'
#' This wrapper function is useful for simulation studies, where the performance
#' of supervised and unsupervised learning methods and algorithms are evaluated.
#' For each data-generating model, we generate \eqn{n_k} observations \eqn{(k =
#' 1, \ldots, K)} from each of \eqn{K} multivariate distributions.
#'
#' Each family returns a list containing a matrix of the multivariate
#' observations generated as well as the class labels for each observation.
#'
#' For details about an individual data-generating family, please see its
#' respective documentation.
#'
#' @param family the family of distributions from which to generate data
#' @param ... optional arguments that are passed to the data-generating function
#' @return named list containing:
#' \describe{
#'   \item{x:}{A matrix whose rows are the observations generated and whose
#'   columns are the \code{p} features (variables)}
#'   \item{y:}{A vector denoting the population from which the observation in
#'   each row was generated.}
#' }
#' @export
#' @examples
#'
#' data_normal <- simdata(family = "normal", n = c(10, 20), mean = c(0, 1), cov = diag(2), seed = 42)
#' data_uniform <- simdata(family = "uniform", delta = 2, seed = 42)
#' data_friedman <- simdata(family = "friedman", experiment = 4, seed = 42)
#' 
simdata <- function(family = c("uniform", "normal", "t", "contaminated", "guo",
                        "friedman"), ...) {
  family <- match.arg(family)
  switch(family,
    normal = simdata_normal(...),
    student = simdata_t(...),
    uniform = simdata_uniform(...),
    contaminated = simdata_contaminated(...),
    guo = simdata_guo(...),
    friedman = simdata_friedman(...)
  )
}
