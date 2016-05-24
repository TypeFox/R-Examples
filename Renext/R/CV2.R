## ****************************************************************************
## DO NOT ROXYGENISE THIS FILE!!! The roxygen doc is was used to
## generate a first draft of Rd files whih have been edited since.
## ****************************************************************************

##' Squared Coefficient of Variation. 
##'
##' Compute the squared Coefficient of Variation of one or several samples provided
##' as a numeric vector or matrix.
##'
##' @title Squared Coefficient of Variation 
##'
##' @param x Numeric vector or matrix. 
##'
##' @return Numeric vector of the squared coefficients of variation.
##' 
##' @note
##' The squared coefficient of variation is the ratio  \eqn{\bar{X}/S^2}{xbar/S^2}
##' where \eqn{\bar{X}}{xbar} and \eqn{S^2} are the sample mean and the sample
##' variance. The variance is computed using the sample size \eqn{n} denominator
##' rather than the usual \eqn{n-1}.
##' 
##' @examples
##' n <- 30; nSamp <- 500
##' X <- matrix(rexp(n * nSamp), nrow= nSamp, ncol = n)
##' W <- CV2(X)
##' plot(density(W), main = "CV2 of exponential samples")
##' 
CV2 <- function(x) {
    
    if (is.null(dim(x))) {
        v <- var(x) * (length(x) - 1) / length(x)
        CV2 <- v / mean(x)^2
       return(CV2)
    } else {
        if (!is.matrix(x)) {
            stop("'x' must be a vector or a matrix")
        }
        n <- ncol(x)
        xbar <- apply(x, 1, mean)
        v <- apply(x, 1, var) * (n - 1) / n
        CV2 <- v / xbar^2
        return(CV2)
    }
    
}

##' Test of exponentiality based on the squared coefficient of
##' variation
##'
##' @title CV2 test of exponentiality
##'
##' @param x Numeric vector giving the sample.
##'
##' @param method Method used to compute the \eqn{p}-value, as
##' in \code{\link{Jackson.test}}.
##'
##' @param nSamp Number of samples used to compute the
##' \eqn{p}-value.
##'
##' @return A list of test results.
##' \item{statistic}{
##' 
##' The test statistic, i.e. the squared coefficient of variation.
##'   
##' }
##' \item{df}{
##' 
##' The sample size.
##' 
##' }
##' \item{p.val}{
##' 
##' The \eqn{p}-value.
##' 
##' }
##' \item{method}{
##' 
##' Description of the test method.
##' 
##' }
##' 
##' @author Yves Deville
##'
##' @note The expectation of CV2 is close to \code{1} for a large sample size 
##' \code{n}.
##' @seealso \code{\link{CV2}}.
##' 
##' @examples
##' n <- 30; nSamp <- 500
##' X <- matrix(rexp(n * nSamp), nrow = nSamp, ncol = n)
##' pVals <- apply(X, 1, function(x) CV2.test(x)$p.value)
##' plot(pVals)  ## should be uniform on (0, 1)
##' 
CV2.test <- function(x,
                     method = c("num", "sim", "asymp"),
                     nSamp = 15000){
    method <- match.arg(method)
    n <- length(x)
    w <- CV2(x)
    if (method == "num") {
        p <- seq(from = 0.01, to = 0.99, by = 0.01)
        table <- qStat(p = p, n = n, type = "Greenwood")
        pVal <- 1 - spline(x = table$q, y = table$p, xout = w)$y
        if (pVal < 0 ) pVal <- 0
        if (pVal > 1 ) pVal <- 1
    } else if (method == "sim") {
        X <- matrix(rexp(n * nSamp), nrow = nSamp)
        W <- CV2(X)
        pVal <- mean(W > w)
    } else if (method == "asymp") {
        wMod <- sqrt(n) * (w - 1) / 2
        pVal <- pnorm(wMod, lower.tail = FALSE)
    }

    pVal <- round(pVal, digits = 3)
    
    list(statistic = w,
         df = n,
         p.value = pVal,
         method = "Wilk's test of exponentiality")
}
