#' Density and log-likelihood of the Bivariate Frank Copula Weibull Count model
#'
#' Compute density and log-likelihood of the Bivariate Frank Copula Weibull
#' Count model.
#'
#' \code{dBivariateWeibullCountFrankCopula} computes the probabilities
#' \eqn{P(X(t) = x(t), Y(t) = y(t))}, where \eqn{X(t),Y(t)} is a bivariate
#' Weibull count process in which the bivariate distribution is modelled by
#' Frank copulas.
#'
#' @param x,y numeric, the desired counts.
#' @param shapeX,shapeY numeric, shape parameters. Either length(x) or
#'     length(1).
#' @param scaleX,scaleY numeric, scale parameters (length(x)).
#' @param theta numeric, Frank copula parameter.
#' @param method character method to be used. Choices are \code{"series_acc"}
#'     (accelerated series expansion) or \code{"conv_dePril"} (convolution by
#'     dePril algorithm).
#' @param time numeric, length of the observation window (defaults to 1).
#' @param conv_steps integer, number of steps to use in the computation of the
#'     integral.
#' @param conv_extrap logical, if \code{TRUE}, Richardson extrapolation will be
#'     applied to improve accuracy.
#' @param series_terms number of terms used in series expansion
#' @param series_acc_niter number of iterations in the acceleration algorithm.
#' @param series_acc_eps double, tolerance to declare convergence in the
#'     acceleration algorithm.
#' @param log TODO
#' @return for \code{dBivariateWeibullCountFrankCopula}, a vector of the
#'     (log-)probabilities
#' @examples
#' ## first 10 cases from "estimationParams.RDS", rounded for presentation
#' gam_weiH <-  0.9530455
#' gam_weiA <-  1.010051
#' theta    <- -0.3703702
#' HG <- c(0, 0, 0, 2, 1, 0, 2, 0, 1, 2)
#' AG <- c(2, 2, 1, 1, 6, 1, 0, 2, 0, 1)
#' lambdaHome <- c(1.5, 1.0, 1.3, 1.8, 1.3, 1.2, 1.3, 1.0, 2.0, 1.4)
#' lambdaAway <- c(1.2, 2.4, 1.3, 0.7, 1.3, 1.4, 0.6, 1.6, 0.6, 1.3)
#'
#' weiFrank0 <- dBivariateWeibullCountFrankCopula(
#'     HG, AG, gam_weiH, lambdaHome, gam_weiA, lambdaAway, theta,
#'     "series_acc", 1, TRUE)
#'
#' weiFrank1 <- dBivariateWeibullCountFrankCopula(
#'     HG, AG, gam_weiH, lambdaHome, gam_weiA, lambdaAway, theta,
#'     "conv_dePril", 1, TRUE, conv_extrap = TRUE)
#'
#'
#' weights <- c(0.01355306, 0.01355306, 0.01355306, 0.01355306, 0.01355306,
#'              0.01355306, 0.01355306, 0.01355306, 0.01357825, 0.01357825)
#'
#' weiFrank2 <- dBivariateWeibullCountFrankCopula_loglik(
#'     HG, AG, gam_weiH, lambdaHome, gam_weiA, lambdaAway, theta,
#'     "conv_dePril", 1, TRUE, conv_extrap = TRUE, weights = weights)
#'
#' weiFrank3 <- dBivariateWeibullCountFrankCopula_loglik(
#'     HG, AG, gam_weiH, lambdaHome, gam_weiA, lambdaAway, theta,
#'     "series_acc", 1, TRUE, weights = weights)
#'
#' cbind(weiFrank0, weiFrank1, weiFrank2, weiFrank3)
#' ## rdname dRenewalFrankCopula_user
#' @export
dBivariateWeibullCountFrankCopula <- function(x, y, shapeX, scaleX,
                                              shapeY, scaleY, theta,
                                              method = c("series_acc",
                                                  "conv_dePril"),
                                              time = 1, log = FALSE,
                                              conv_steps = 100,
                                              conv_extrap = TRUE,
                                              series_terms = 50,
                                              series_acc_niter = 300,
                                              series_acc_eps = 1e-10) {

    method <- match.arg(method)
    n <- length(x)
    stopifnot(length(y) == n,
              length(scaleX) == n,
              length(scaleY) == n)

    uni <- FALSE
    if (length(shapeX) == 1 & length(shapeY) == 1)
        uni <- TRUE
    else if (length(shapeX) == n & length(shapeY) == n)
        uni <- FALSE
    else
        stop("check shapeX and shapeY: should be length(1) or length(x) !")

    switch(method,
           "series_acc" = {
               if (uni)
                   return(
                       dWeibullInterArrivalCountFrankCopula_uni(
                           x, y, shapeX, shapeY, scaleX, scaleY, theta,
                           time, log, series_terms, series_acc_niter,
                           series_acc_eps)
                       )
               else
                  return(
                       dWeibullInterArrivalCountFrankCopula(
                           x, y, shapeX, shapeY, scaleX, scaleY, theta,
                           time, log, series_terms, series_acc_niter,
                           series_acc_eps)
                       )
           },
           "conv_dePril" = {
               if (length(shapeX) == 1)
                   shapeX <- rep(shapeX, n)
               if (length(shapeY) == 1)
                   shapeY <- rep(shapeY, n)

               .createDistPar <- function(ind, Xflag = TRUE) {
                   if (Xflag)
                       list(shape = as.numeric(shapeX[ind]),
                            scale = as.numeric(scaleX[ind])
                            )
                   else
                       list(shape = as.numeric(shapeY[ind]),
                            scale = as.numeric(scaleY[ind])
                            )
               }
               distPX <- lapply(seq(along = x),  .createDistPar ,
                                Xflag = TRUE)
               distPY <- lapply(seq(along = x),  .createDistPar ,
                                Xflag = FALSE)

               return(dRenewalFrankCopula_bi(x, y, "weibull", "weibull",
                                             distPX, distPY,
                                             theta, time, log,
                                             conv_steps, conv_extrap)
                      )
           })
}

#' % Log-likelihood of the Bivariate Frank Copula Weibull Count model
#'
#' @param na.rm logical, should \code{NA}s (obtained from log of small
#'     probabilities) be replaced with the smallest allowed probability?
#' @param weights numeric vector of weights to apply. If \code{NULL}, a vector
#'     of ones.
#' @return for \code{dBivariateWeibullCountFrankCopula_loglik}, the
#'     log-likelihood of the model, a number
#' @rdname dBivariateWeibullCountFrankCopula
#' @export
dBivariateWeibullCountFrankCopula_loglik <- function(x, y, shapeX, scaleX,
                                                     shapeY, scaleY, theta,
                                                     method = c("series_acc",
                                                         "conv_dePril"),
                                                     time = 1, na.rm = TRUE,
                                                     conv_steps = 100,
                                                     conv_extrap = TRUE,
                                                     series_terms = 50,
                                                     series_acc_niter = 300,
                                                     series_acc_eps = 1e-10,
                                                     weights = NULL) {
    method <- match.arg(method)
    if (is.null(weights))
        weights <- rep(1, length(x))

    pbs <- dBivariateWeibullCountFrankCopula(x, y, shapeX, scaleX,
                                             shapeY, scaleY,
                                             theta, method, time, TRUE,
                                             conv_steps, conv_extrap,
                                             series_terms, series_acc_niter,
                                             series_acc_eps) * weights
    if (na.rm)
        pbs[is.na(pbs)] <- .logNaReplace()

    sum(pbs)
}
