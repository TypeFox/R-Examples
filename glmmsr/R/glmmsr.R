#' glmmsr: fit GLMMs with various approximation methods
#'
#' The glmmsr package provides functions to conduct inference about generalized linear
#' mixed models, giving the user a choice about which method to use to
#' approximate the likelihood.
#'
#' In addition to the Laplace and adaptive Gaussian quadrature approximations,
#' which are borrowed from \code{lme4}, the likelihood may
#' be approximated by the sequential reduction approximation
#' or an importance sampling approximation. These methods
#' provide an accurate approximation to the likelihood in some situations
#' where it is not possible to use adaptive Gaussian quadrature.
#'
#' The main function of the glmmsr package is \code{\link{glmm}}, which is
#' used to fit the GLMM. Its interface allows a larger class of models than
#' those allowed by \code{lme4}, including
#' structured pairwise comparison models.
#'
#' @references Helen E. Ogden (2015). A sequential reduction method for
#'   inference in generalized linear mixed models. Electronic Journal of
#'   Statistics 9: 135-152. doi:
#'   \href{http://dx.doi.org/10.1214/15-EJS991}{10.1214/15-EJS991}
#' @docType package
#' @name glmmsr
#' @useDynLib glmmsr
#' @importFrom Rcpp sourceCpp
#' @importFrom methods as
#' @importFrom stats dnorm drop.terms family formula gaussian glm lm model.matrix
#'  model.response model.weights na.fail offset optim pnorm rnorm terms update.formula
NULL

#' A dataset simulated from a two-level model
#'
#' @example inst/examples/two_level.R
"two_level"

#' A dataset simulated from a three-level model
#'
#' @example inst/examples/three_level.R
"three_level"

.onUnload <- function (libpath) {
  library.dynam.unload("glmmsr", libpath)
}
