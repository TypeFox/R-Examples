##' Plotting positions for exponential retun level plots. 
##'
##' The plotting positions are numeric values to use as the abscissae
##' corresponding to the order statistics in an exponential return
##' level plot. They range from 1 to about \eqn{\log n}{log(n)}.
##' They can be related to the plotting positions given by
##' \code{\link{ppoints}}.
##' 
##' @title Plotting positions for exponential return levels
##' @param n sample size.
##' @return numeric vector of plotting positions with length
##' \code{n}.
##' @author Yves Deville
##' @seealso \code{\link{ppoints}}.
##' @details
##' The computed values  
##' deqn{H_{i} = \frac{1}{n} + \frac{1}{n-1} + \dots + \frac{1}{n + 1 -i}}{
##' H[i] = 1 / n + 1 / (n + 1) + ... + 1 / (n + 1 - i)}
Hpoints <- function(n) {
    if (length(n) > 1L) n <- length(n)
    if (n == 0) return(numeric(0))
    cumsum(1 / rev(seq(1L:n)))
}
