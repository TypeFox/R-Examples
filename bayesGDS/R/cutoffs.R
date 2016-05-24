#' @name get.cutoffs
#' @aliases get.cutoffs
#' @title Draw thresholds for the accept-reject stage of the BD
#' sampling algorithm.
#' @description Returns a vector of log(u), where u is the threshold
#' to determine if a proposal draw should be accepted as a draw from
#' the target posterior distribution.
#' @param log.phi Vector of log.phi from the proposal draws.  All must
#' be non-positive.
#' @param n.draws an integer.  number of draws to be taken from the
#' target posterior.
#' @return a numeric vector for v = -log.u (the thresholds for the
#' accept-reject stage).
#' @details For use in conjunction with the Braun and Damien (2012)
#' Generalized Direct Sampling algorithm.  This is usually not called
#' directly (and, thus, it is not exported), since it is called from
#' the sample.GDS function.
get.cutoffs <- function(log.phi, n.draws) {

    if (any(log.phi>0) ) stop("all values of log phi must be non-positive")

    M <- length(log.phi)
    v <- c(-log.phi[order(log.phi,decreasing=TRUE)],Inf)
    v1 <- min(v)
    p.int <- as.numeric(exp(log(seq(1,M)) -v1
                            + log(exp(v1-v[1:M]) - exp(v1-v[2:(M+1)]))))
    j <- sample(1:M,n.draws,replace=TRUE,prob=p.int)
    u <- runif(n.draws)
    d <- as.numeric(v[j] - log(1-u + u*exp(-(v[j+1]-v[j]))))


    return(d)
}

