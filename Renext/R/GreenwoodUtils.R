## ****************************************************************************
## DO NOT ROXYGENISE THIS FILE!!! The roxygen doc is was used to
## generate a first draft of Rd files whih have been edited since.
## ****************************************************************************

##' Probability that the Greenwood's statistic is smaller than one.
##'
##' The probability was computed by using the approximation of the
##' quantile function of the GreenWood's statistic returned by
##' \code{\link{qStat}}. The result is found by interpolating the
##' distribution function for \eqn{x = 1}.
##'
##' @title Probability that the Greenwood's statistic is smaller than one
##'
##' @param n Sample size.
##'
##' @return Probability that the Greenwood's statistic is smaller than
##' one. For a random sample of an exponential distribution with size
##' \eqn{n}, this is the probability that the coefficient of variation
##' is less than one, or the probability that the ML estimate of the
##' GPD shape parameter \eqn{\xi} is negative.
##' 
##' @author Yves Deville
##'
##' @examples
##' n <- 8:500
##' plot(n, pGreenwood1(n), type = "l", col = "orangered", lwd = 2,
##'      log ="x", ylim =c(0.5, 0.7), main = "slow convergence to 0.5")
##' grid() ; abline(h = 0.5, col = "SpringGreen")
pGreenwood1 <- function(n) {
    if (any(n < min(pGreenwoodData$n)) || any(n > max(pGreenwoodData$n))){
        stop("'n' must be in the range from ", min(pGreenwoodData$n),
             " to ", max(pGreenwoodData$n))
    }
    spline(x = pGreenwoodData$n, y = pGreenwoodData$prob, xout = n)$y
}

##' Quantile of a test statistic.
##'
##' @details
##' The function provides an approximation of the distribution for several
##' statistics.
##' 
##' \itemize{
##' 
##'  \item For \code{"Greenwood"}, the statistic is Greenwood's one. The
##' distribution is that of the squared coefficient of variation of a
##' sample of size \code{n} from the exponential distribution.
##'
##' \item For \code{"Jackson"}, the statistic is Jackson's satistic,
##' see \code{\link{Jackson}}.
##'
##' \item For \code{"logLRGPD"} and \code{"logLRLomax"}, the statistic 
##' is the log of the likelihood ratio of a sample from the exponential
##' distribution. The log-likelihoods are for an exponential
##' distribution compared to a GPD with non-zero shape, or to a GPD
##' with \emph{positive shape} (equivalently, a Lomax distribution).
##' 
##' \item For \code{"logLRGEV"} and \code{"logLRFrechet"}, the statistic 
##' is the log of the likelihood ratio of a sample from the Gumbel
##' distribution. The log-likelihoods are for a Gumbel
##' distribution compared to a GEV with non-zero shape, or to a GEV
##' with \emph{positive shape} (equivalently, a Frechet distribution).
##' 
##' }
##' 
##' For log of likelihood ratios, these are multiplied by \code{2}, so that
##' they compare to a chi-square statistic with one degree of freedom.
##' 
##' @title Quantiles of a test statistic 
##'
##' @param p Numeric vector of probabilities. Very small values
##' (\code{p < 0.01}) or very large ones (\code{p > 0.99}) will be
##' truncated to maintain a realistic level of precision.
##'
##' @param n Sample size.
##'
##' @param type The type of satistic, see \bold{Details}.
##'
##' @param outNorm Logical. If \code{TRUE} the output is normalized in
##' a such fashion that its distribution is the asympotic one (i.e.
##' standard normal in practice). When \code{FALSE}, the quantiles
##' are given in the true scale of the statisic: \eqn{\textrm{CV}^2}{CV^2}, Jackson.
##' For LR statistics this argument has no impact.
##'
##' @return A vector of quantiles.
##'
##' @author Yves Deville
##'
##' @note The precision of the result given is limited, and is about
##' two-digits. This function is not intented to be used as such and is only
##' provided for information.
##' 
##' @examples
##' res <- qStat(n = 40, type = "Greenwood")
##' plot(res$q, res$p, type = "o")
qStat <- function(p, n,
                  type = c("Greenwood", "Jackson", "logLRGPD", "logLRLomax",
                      "logLRGEV", "logLRFrechet"),
                  outNorm = FALSE) {

    type <- match.arg(type)
    object <- quantApp[[type]]
    
    if (missing(p)) p <- object$pGrid
    
    if (length(n) > 1L) stop("'n' must be of length 1")
    if (n < object$n0) stop(sprintf("n must be >= %d", object$n0))

    ind <- (p < 0.005)
    if (any(ind)) {
        warning("'p' contains values smaller than 0.005. Truncated.")
        p[ind] <- 0.005
    }
    ind <- (p > 0.995)
    if (any(ind)) {
        warning("'p' contains values larger than 0.995. Truncated.")
        p[ind] <- 0.995
    }
        
    q <- qnorm(p)
    X <- splineDesign(knots = object$qKnots, ord = object$ord, x = q)
    
    if (object$type == "logLRLomax") {
        X[p <= pGreenwood1(n), ] <- 0
    }

    if (n < object$n1) {
        alpha <- object$Alpha[object$nA ==n, ]
        res <- X %*% alpha
    } else {
        b <- t(object$Beta) %*% n^(c(0, -1, -2, -3) / 2)
        res <-  X %*% b
    }

    if (!outNorm && object$normalize) {
        res <- object$gamma + object$delta * res / sqrt(n) 
    }
    
    list(p = p, q = as.numeric(res))
}
