#' Compute count probabilities using convolution (bi)
#'
#' Compute count probabilities using one of several convolution methods for the
#' distributions with builtin support in this package.
#'
#' The following convolution methods are implemented: "dePril", "direct", and
#' "naive".
#'
#' The builtin distributions currently are Weibull, gamma, generalised gamma and
#' Burr.
#'
#' @param x integer (vector), the desired count values.
#' @inheritParams surv
#' @param nsteps unsiged integer number of steps used to compute the integral.
#' @param time double time at wich to compute the probabilities. Set to 1 by
#' default.
#' @param extrap logical if \code{TRUE}, Richardson extrapolation will be
#' applied to improve accuracy.
#' @param log logical if \code{TRUE} the log-probability will be returned.
#' @param method TODO
#' @return vector of probabilities P(x(i)) for i = 1, ..., n where n is
#' \code{length} of \code{x}.
#' @examples
#' x <- 0:10
#' lambda <- 2.56
#' p0 <- dpois(x, lambda)
#' ll <- sum(dpois(x, lambda, TRUE))
#'
#' err <- 1e-6
#' ## all-probs convolution approach
#' distPars <- list(scale = lambda, shape = 1)
#' pmat_bi <- dCount_conv_bi(x, distPars, "weibull", "direct",
#'                           nsteps = 200)
#'
#' ## user pwei
#' pwei_user <- function(tt, distP) {
#'     alpha <- exp(-log(distP[["scale"]]) / distP[["shape"]])
#'     pweibull(q = tt, scale = alpha, shape = distP[["shape"]],
#'              lower.tail = FALSE)
#' }
#'
#' pmat_user <- dCount_conv_user(x, distPars, c(1, 2), pwei_user, "direct",
#'                               nsteps = 200)
#' max((pmat_bi - p0)^2 / p0)
#' max((pmat_user - p0)^2 / p0)
#'
#' ## naive convolution approach
#' pmat_bi <- dCount_conv_bi(x, distPars, "weibull", "naive",
#'                           nsteps = 200)
#' pmat_user <- dCount_conv_user(x, distPars, c(1, 2), pwei_user, "naive",
#'                               nsteps = 200)
#' max((pmat_bi- p0)^2 / p0)
#' max((pmat_user- p0)^2 / p0)
#'
#' ## dePril conv approach
#' pmat_bi <- dCount_conv_bi(x, distPars, "weibull", "dePril",
#'                           nsteps = 200)
#' pmat_user <- dCount_conv_user(x, distPars, c(1, 2), pwei_user, "dePril",
#'                               nsteps = 200)
#' max((pmat_bi- p0)^2 / p0)
#' max((pmat_user- p0)^2 / p0)
#'
#' @export
dCount_conv_bi <- function(x, distPars,
                           dist = c("weibull", "gamma", "gengamma",
                               "burr"),
                           method = c( "dePril", "direct", "naive"),
                           nsteps = 100, time = 1.0,
                           extrap = TRUE, log = FALSE) {

    dist <- match.arg(dist)
    method <- match.arg(method)

    switch(method,
           "direct" = {
               return(dCount_allProbs_bi(x, distPars, dist,
                                         nsteps, time, extrap, log))
           },
           "naive" = {
               return(dCount_naive_bi(x, distPars, dist,
                                      nsteps, time, extrap, 0, log))
           },
           "dePril" = {
               return(dCount_dePril_bi(x, distPars, dist,
                                       nsteps, time, extrap, 0, log))
           }
           )
}

#' Compute count probabilities using convolution (user)
#'
#' Compute count probabilities using one of the convolution methods for user
#' defined survival functions.
#'
#' @param extrapolPars ma::vec of length 2. The extrapolation values.
#' @param survR Rcpp::Function user passed survival function; should have the
#' signature \code{function(t, distPars)} where \code{t} is a real number (>0)
#' where the survival function is evaluated and \code{distPars} is a list of
#' distribution parameters. It should return a double value.
#' @inheritParams dCount_conv_bi
#' @return vector of probabilities P(x(i)) for i = 1, ..., n where n is
#' \code{length} of \code{x}.
#' @examples
#' ## see examples for dCount_conv_bi
#' @export
dCount_conv_user <- function(x, distPars, extrapolPars, survR,
                             method = c( "dePril", "direct", "naive"),
                             nsteps = 100, time = 1.0,
                             extrap = TRUE, log = FALSE) {

    method <- match.arg(method)
    switch(method,
           "direct" = {
               return(dCount_allProbs_user(x, distPars, extrapolPars, survR,
                                           nsteps, time, extrap, log))
           },
           "naive" = {
               return(dCount_naive_user(x, distPars, extrapolPars, survR,
                                        nsteps, time, extrap, 0, log))
           },
           "dePril" = {
               return(dCount_dePril_user(x, distPars, extrapolPars, survR,
                                         nsteps, time, extrap, 0, log))
           }
           )
}

