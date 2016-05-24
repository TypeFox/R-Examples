#' wrapper to the weibull count probability
#'
#' wrapper to the univariate weibull count probability
#'
#' User could choose between one of the following methods:
#' \describe{
#' \item{series_mat}{series expansion using matrix techniques.}
#' \item{series_acc}{Euler-van Wijngaarden accelerated series expansion}
#' \item{conv_direct}{direct convolution method of section 2.}
#' \item{conv_naive}{naive convolurion described in section 3.1}
#' \item{conv_dePril}{dePril convolution described in section 3.2}
#' }
#' The valiue of the inputs can be left to their default value if user
#' is not aware of the algorithm details.
#' @param scale numeric (length 1), scale parameter of the Weibull count.
#' @param shape numeric (length 1), shape parameter of the Weibull count.
#' @param x integer (vector), the desired count values.
#' @param time double, length of the observation window (defaults to 1).
#' @param log logical, if TRUE, the log of the probability will be returned.
#' @param method character one of the available methods. See details.
#' @param conv_steps numeric number of steps used for the extrapolation.
#' @param conv_extrap logical should Richardson extrappolation be applied ?
#' @param series_terms numeric number of terms in the series expansion.
#' @param series_acc_niter numeric number of iteration in the
#' Euler-van Wijngaarden algorithm.
#' @param series_acc_eps numeric tolerance of convergence in the
#' Euler-van Wijngaarden algorithm.
#' @return vector of probabilities \eqn{P(x(i)), i = 1, \dots n} where
#' \code{n = length(x)}.
#' @export
dWeibullCount <- function(x, shape, scale,
                          method = c("series_acc", "series_mat",
                              "conv_direct", "conv_naive", "conv_dePril"),
                          time = 1, log = FALSE, conv_steps = 100,
                          conv_extrap = TRUE, series_terms = 50,
                          series_acc_niter = 300, series_acc_eps = 1e-10) {

    method <- match.arg(method)
    distPars <- list(scale = scale, shape = shape)
    dist <- "weibull"
    
    switch(method,
           "conv_direct" = {
               return(dCount_allProbs_bi(x, distPars, dist,
                                         conv_steps, time, conv_extrap, log)
                      )
           },
           "conv_naive" = {
               return(dCount_naive_bi(x, distPars, dist,
                                      conv_steps, time, conv_extrap, 0, log))
           },
           "conv_dePril" = {
               return(dCount_dePril_bi(x, distPars, dist,
                                       conv_steps, time, conv_extrap, 0, log))
           },
           "series_mat" = {
               return(dWeibullCount_mat(x, shape, scale, time, log,
                                        series_terms))
           },
           "series_acc" = {
               return(dWeibullCount_acc(x, shape, scale, time, log,
                                        series_terms, series_acc_niter,
                                        series_acc_eps)
                   )
           }
           )
}


#' wrapper to the weibull count log-likelihood
#'
#' wrapper to the univariate weibull log-likelihood
#'
#' @param na.rm logical if TRUE, the \code{NA} (produced by taking the log of
#' very small probabilities) will be replaced by the smallest allowed probaility;
#' default = \code{TRUE}
#' @param weights numeric vector of weights to apply. If \code{NULL}, one will
#' be applied.
#' @inheritParams dWeibullCount
#' @return double log-likelihood of the count process
#' @rdname dWeibullCount
#' @export
dWeibullCount_loglik <- function(x, shape, scale, 
                                 method = c("series_acc", "series_mat",
                                     "conv_direct", "conv_naive", "conv_dePril"),
                                 time = 1, na.rm = TRUE, conv_steps = 100,
                                 conv_extrap = TRUE, series_terms = 50,
                                 series_acc_niter = 300,
                                 series_acc_eps = 1e-10, weights = NULL) {
    if (is.null(weights))
        weights <- rep(1, length(x))
    
    method <- match.arg(method)
    
    ## -------------- rescale parameters if needed
    if (length(x) != length(scale))
        stop("length x should be the same as length scale !")

    if (length(shape) == 1)
        shape <- rep(shape, length(x))
    else if (length(x) != length(shape))
        stop("shape should be length 1 or same as length x !")

    .getDistPars <- function(ind)
        list(scale = scale[ind], shape = shape[ind])
    
    distPars <- lapply(seq(along = x), .getDistPars)
    dist <- "weibull"
    
    switch(method,
           "conv_direct" = {
               .computeOneProb <- function(ind)
                   dCount_allProbs_scalar_bi(x[ind], distPars[[ind]],
                                             dist, conv_steps, time,
                                             conv_extrap, logFlag = TRUE)
           },
           "conv_naive" = {
               .computeOneProb <- function(ind)
                   dCount_naive_scalar_bi(x[ind], distPars[[ind]],
                                          dist, conv_steps, time,
                                          conv_extrap, logFlag = TRUE)

           },
           "conv_dePril" = {
               .computeOneProb <- function(ind)
                   dCount_dePril_scalar_bi(x[ind], distPars[[ind]],
                                           dist, conv_steps, time,
                                           conv_extrap, logFlag = TRUE)   
           },
           "series_mat" = {
               .computeOneProb <- function(ind)
                   dWeibullCount_mat(x[ind], shape[ind], scale[ind],
                                     time, TRUE,
                                     series_terms)
           },
           "series_acc" = {
               .computeOneProb <- function(ind)
                   dWeibullCount_acc(x[ind], shape[ind], scale[ind],
                                     time, TRUE,
                                     series_terms, series_acc_niter,
                                     series_acc_eps)
               
           })

     pbs <- sapply(seq(along = x), .computeOneProb) * weights
     if (na.rm)
         pbs[is.na(pbs)] <- .logNaReplace()
     
     sum(pbs)
}


#' Weibull count expected value and variance
#'
#' Expcted value and varinace of the weibull count process
#'
#' @param xmax unsigned integer maximum count to be used.
#' @inheritParams dWeibullCount
#' @return double log-likelihood of the count process
#' @rdname dWeibullCount
#' @export
evWeibullCount <- function(xmax, shape, scale, 
                           method = c("series_acc", "series_mat",
                               "conv_direct", "conv_naive", "conv_dePril"),
                           time = 1, conv_steps = 100,
                           conv_extrap = TRUE, series_terms = 50,
                           series_acc_niter = 300,
                           series_acc_eps = 1e-10) {

    method <- match.arg(method)

    x <- 1:xmax
    px <- dWeibullCount(x, shape, scale, method, time, FALSE,
                        conv_steps, conv_extrap, series_terms,
                        series_acc_niter, series_acc_eps)
    
    ev <- sum(x * px)
    ev2 <- sum(x^2 * px)
    var <- ev2 - ev^2
    list(ExpectedValue = ev, Variance = var)
}
