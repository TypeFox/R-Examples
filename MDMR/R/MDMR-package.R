#' Multivariate Distance Matrix Regression
#'
#' \code{MDMR} allows a user to conduct multivariate distance matrix regression
#' using analytic p-values and measures of effect size described by McArtor &
#' Lubke (submitted). Analytic p-values are computed using the R package
#' CompQuadForm (Duchesne & De Micheaux, 2010).
#'
#' @section Usage:
#'  To access this package's tutorial, type the following line into the console:
#'
#'  \code{vignette('mdmr-vignette')}
#'
#'  There are two primary functions that comprise this package:
#'  \code{\link{mdmr}}, which regresses a distance matrix onto a set of
#'  predictors, and \code{\link{delta}}, which computes measures of univariate
#'  effect size in the context of multivariate distance matrix regression. The
#'  help files of both functions provide more general information than the
#'  package vignette.
#'
#' @references Davies, R. B. (1980). The Distribution of a Linear Combination of
#'  chi-square Random Variables. Journal of the Royal Statistical Society.
#'  Series C (Applied Statistics), 29(3), 323-333.
#'
#'  Duchesne, P., & De Micheaux, P.L. (2010). Computing the distribution of
#'  quadratic forms: Further comparisons between the Liu-Tang-Zhang
#'  approximation and exact methods. Computational Statistics and Data
#'  Analysis, 54(4), 858-862.
#'
#'  McArtor, D.B. & Lubke, G.H. (submitted). Extending multivariate distance
#'  matrix regression with an effect size measure and the distribution of the
#'  test statistic.
#'
#' @examples
#'data(mdmrdata)
#'D <- dist(Y.mdmr, method = 'euclidean')
#'
#'mdmr.res <- mdmr(X = X.mdmr, D = D)
#'summary(mdmr.res)
#'
#'mdmr.delta <- delta(X = X.mdmr, Y = Y.mdmr, dtype = 'euclidean',
#'                    niter = 1, seed = 12345)
#'
#' @importFrom CompQuadForm davies
#' @importFrom  parallel mclapply
#'
#' @docType package
#' @name MDMR-package
#' @aliases MDMR
NULL
