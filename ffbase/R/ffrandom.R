
#' Generate \code{ff} vector with draws from distribution
#'
#' A convenience function to generate \code{ff} vectors with draws from random
#' distributions using functions such as \code{\link{runif}},
#' \code{\link{rnorm}} and \code{\link{rlnorm}}. 
#'
#' @param n number of observations
#' @param rfun a function generating the draws from the random distribution. 
#'   This function should expect the number of draws as its first argument. 
#'   Valid examples are the routines \code{\link{runif}}, \code{\link{rnorm}},
#'   and \code{\link{rlnorm}}. 
#' @param ... additional arguments are passed on to \code{rfun}.
#' @param vmode the vmode of the resulting vector. See \code{\link{ff}}. If 
#'   none given the vmode is determined from a single draw from \code{rfun}. 
#'
#' @details Before generating the vector a single draw is taken from the 
#' distribution. This might be important if one tries to reproduce draws 
#' directly from \code{rfun}. 
#' 
#' @return An \code{ff} vector with draws from the distribution.
#'
#' @example ../examples/ffrandom.R
#'
#' @export
#'

ffrandom <- function(n, rfun = runif, ..., vmode = NULL) {
    r <- ff(rfun(1), length=n, vmode=vmode)
    for (i in chunk(r)) {
        Log$chunk(i)
        ni <- diff(range(i))+1
        r[i] <- rfun(ni, ...)
    }
    r
}

