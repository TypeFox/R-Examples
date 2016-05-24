#' Wrapper function to generate data from a variety of data-generating models.
#'
#' We provide a wrapper function to generate from three data-generating models:
#' \describe{
#'   \item{\code{\link{sim_unif}}}{Five multivariate uniform distributions}
#'   \item{\code{\link{sim_normal}}}{Multivariate normal distributions with intraclass covariance matrices}
#'   \item{\code{\link{sim_student}}}{Multivariate Student's t distributions each with a common covariance matrix}
#' }
#'
#' For each data-generating model, we generate \eqn{n_m} observations
#' \eqn{(m = 1, \ldots, M)} from each of \eqn{M} multivariate distributions so
#' that the Euclidean distance between each of the population centroids and the
#' origin is equal and scaled by \eqn{\Delta \ge 0}. For each model, the argument
#' \code{delta} controls this separation.
#'
#' This wrapper function is useful for simulation studies, where the efficacy of
#' supervised and unsupervised learning methods and algorithms are evaluated as a
#' the population separation is increased.
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
#' set.seed(42)
#' uniform_data <- sim_data(family = "uniform")
#' normal_data <- sim_data(family = "normal", delta = 2)
#' student_data <- sim_data(family = "student", delta = 1, df = 1:5)
sim_data <- function(family = c("uniform", "normal", "student"), ...) {
  family <- match.arg(family)
  switch(family,
    uniform = sim_unif(...),
    normal = sim_normal(...),
    student = sim_student(...)
  )
}
