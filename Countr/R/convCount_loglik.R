#' Log-likelihood of a count probability computed by convolution (bi)
#'
#' Compute the log-likelihood of a count model using  convolution
#' methods to compute the probabilities.
#' \code{dCount_conv_loglik_bi} is for the builtin distributions.
#' \code{dCount_conv_loglik_user} is for user defined survival functions.
#'
#' @param distPars list of the same length as x with each slot being itself a
#'     named list containing the distribution parameters corresponding to
#'     \code{x[i]}
#' @param method character, convolution method to be used; choices are
#'     \code{"dePril"} (section 3.2), \code{"direct"} (section 2) or
#'     \code{"naive"} (section 3.1).
#' @inheritParams dCount_allProbs_bi
#' @param na.rm logical, if TRUE, the \code{NA}s (produced by taking the log of
#'     very small probabilities) will be replaced by the smallest allowed
#'     probability; default is \code{TRUE}.
#' @param weights numeric vector of weights to apply. If \code{NULL}, a vector
#'     of ones.
#' @return numeric, the log-likelihood of the count process
#' @examples
#' x <- 0:10
#' lambda <- 2.56
#' distPars <- list(scale = lambda, shape = 1)
#' distParsList <- lapply(seq(along = x), function(ind) distPars)
#' extrapolParsList <- lapply(seq(along = x), function(ind) c(2, 1))
#' ## user pwei
#' pwei_user <- function(tt, distP) {
#'     alpha <- exp(-log(distP[["scale"]]) / distP[["shape"]])
#'     pweibull(q = tt, scale = alpha, shape = distP[["shape"]],
#'              lower.tail = FALSE)
#' }
#'
#' ## log-likehood allProbs Poisson
#' dCount_conv_loglik_bi(x, distParsList,
#'                       "weibull", "direct", nsteps = 400)
#'
#' dCount_conv_loglik_user(x, distParsList, extrapolParsList,
#'                         pwei_user, "direct", nsteps = 400)
#'
#' ## log-likehood naive Poisson
#' dCount_conv_loglik_bi(x, distParsList,
#'                       "weibull", "naive", nsteps = 400)
#'
#' dCount_conv_loglik_user(x, distParsList, extrapolParsList,
#'                         pwei_user, "naive", nsteps = 400)
#'
#' ## log-likehood dePril Poisson
#' dCount_conv_loglik_bi(x, distParsList,
#'                       "weibull", "dePril", nsteps = 400)
#'
#' dCount_conv_loglik_user(x, distParsList, extrapolParsList,
#'                         pwei_user, "dePril", nsteps = 400)
#' @export
dCount_conv_loglik_bi <- function(x, distPars,
                                  dist = c("weibull", "gamma", "gengamma",
                                      "burr"),
                                  method = c( "dePril", "direct", "naive"),
                                  nsteps = 100,
                                  time = 1.0, extrap = TRUE,
                                  na.rm = TRUE, weights = NULL) {
    dist <- match.arg(dist)
    method <- match.arg(method)
    if (is.null(weights))
        weights <- rep(1, length(x))

    ## case length(x) == 1 treated separatly
    if (length(x) == 1) {
        ## distPars
        if ((length(distPars) > 1) | !is.null(names(distPars))) {
            temp <- list()
            temp[[1]] <- distPars
            distPars <- temp
        }
    }

    ## check inputs length
    if (length(x) != length(distPars))
        stop("length x should be the same as length distPars !")

    switch(method,
           "direct" = {
               .computeOneProb <- function(ind)
                   dCount_allProbs_scalar_bi(x[ind], distPars[[ind]],
                                             dist, nsteps, time,
                                             extrap, logFlag = TRUE)
           },
           "naive" = {
               .computeOneProb <- function(ind)
                   dCount_naive_scalar_bi(x[ind], distPars[[ind]],
                                          dist, nsteps, time,
                                          extrap, logFlag = TRUE)

           },
           "dePril" = {
               .computeOneProb <- function(ind)
                   dCount_dePril_scalar_bi(x[ind], distPars[[ind]],
                                           dist, nsteps, time,
                                           extrap, logFlag = TRUE)
           }
           )

    pbs <- sapply(seq(along = x), .computeOneProb) * weights
    if (na.rm)
        pbs[is.na(pbs)] <- .logNaReplace()

    sum(pbs)
}

## Log-likelihood of a count probability computed by convolution (user)
##log-likelihood of a count probability computed by convolution methods for user
##passed survival functions.
#' @param extrapolPars list of same length as x where each slot is a vector of
#'     length 2 (the extrapolation values to be used) corresponding to
#'     \code{x[i]}.
#' @param survR a user defined survival function; should have the signature
#'     \code{function(t, distPars)} where \code{t} is a real number (>0) where
#'     the survival function is evaluated and \code{distPars} is a list of
#'     distribution parameters. It should return a double value.
#' @inheritParams dCount_conv_loglik_bi
## @return double, the log-likelihood of the count process
#' @rdname dCount_conv_loglik_bi
#' @examples
#' ## see dCount_conv_loglik_bi()
#' @export
dCount_conv_loglik_user <- function(x, distPars, extrapolPars, survR,
                                    method = c( "dePril", "direct", "naive"),
                                    nsteps = 100, time = 1.0, extrap = TRUE,
                                    na.rm = TRUE, weights = NULL) {
     if (is.null(weights))
         weights <- rep(1, length(x))

    ## case length(x) == 1 treated separatly
    if (length(x) == 1) {
        ## distPars
        if ((length(distPars) > 1) | !is.null(names(distPars))) {
            temp <- list()
            temp[[1]] <- distPars
            distPars <- temp
        }
    }

    ## check inputs length
    if (length(x) != length(distPars))
        stop("length x should be the same as length distPars !")
    if (length(x) != length(extrapolPars))
        stop("length x should be the same as length extrapolPars !")

    switch(method,
           "direct" = {
               .computeOneProb <- function(ind)
                   dCount_allProbs_scalar_user(x[ind], distPars[[ind]],
                                               extrapolPars[[ind]],
                                               survR, nsteps, time,
                                               extrap, logFlag = TRUE)
           },
           "naive" = {
               .computeOneProb <- function(ind)
                   dCount_naive_scalar_user(x[ind], distPars[[ind]],
                                            extrapolPars[[ind]],
                                            survR, nsteps, time,
                                            extrap, logFlag = TRUE)

           },
           "dePril" = {
               .computeOneProb <- function(ind)
                   dCount_dePril_scalar_user(x[ind], distPars[[ind]],
                                             extrapolPars[[ind]],
                                             survR, nsteps, time,
                                             extrap, logFlag = TRUE)
           }
           )

     pbs <- sapply(seq(along = x), .computeOneProb) * weights
     if (na.rm)
         pbs[is.na(pbs)] <- .logNaReplace()

     sum(pbs)
}

