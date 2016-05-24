## ****************************************************************************
## DO NOT ROXYGENISE THIS FILE!!! The roxygen doc is was used to
## generate a first draft of Rd files whih have been edited since.
## ****************************************************************************


##' Likelihood-Ratio statistic for the exponential vs. GPD.
##'
##' The Likelihood-Ratio statisitc is actually \eqn{W:=-2 \log
##' \textrm{LR}}{-2 log LR} where LR is the ratio of the likelihoods
##' \emph{exponential} to \emph{alternative distribution}. 
##' @title  Likelihood-Ratio statistic
##'
##' @param x Numeric vector of sample values.
##' 
##' @param alternative Character string describing the alternative hypothesis
##' 
##' @return The LR statistic value.
##'
##' @note When the alternative is \code{"lomax"} or \code{"maxlo"}
##' the statistic has a distribution of mixed type under the null
##' hypothesis of exponentiality. This is a mixture of a distribution
##' of continuous type and of a Dirac mass at LR = 0.
LRExp <- function(x,
                  alternative = c("lomax", "GPD", "gpd", "maxlo")){
    
    alternative <- match.arg(alternative)
    n <- length(x)
    fit <- fGPD(x, cov = FALSE)
    if (alternative == "lomax") {
        if (fit$CV <= 1.0) {
            W <- 0
        } else {
            W <- 2 * ( fit$loglik + n * log(mean(x)) + n )
        }
    } else if (alternative %in% c("GPD", "gpd")) {
        W <- 2 * ( fit$loglik + n * log(mean(x)) + n )
    } else if (alternative == "maxlo") {
        if (fit$CV >= 1.0) {
            W <- 0
        } else {
            W <- 2 * ( fit$loglik + n * log(mean(x)) + n )
        }
    }
    W
}
    
##' Likelihood Ratio test of exponentiality vs. GPD.
##'
##' The \emph{asymptotic} distribution of the Likelihood-ratio
##' statistic is known. For the GPD alternative, this is a chi-square
##' distribution with one df.  For the Lomax alternative, this is the
##' distribution of a product \eqn{XY} where \eqn{X} and \eqn{Y} are
##' two independent random variables following a Bernouilli
##' distribution with probability parameter \eqn{p = 0.5} and a
##' chi-square distribution with one df.
##'
##' When \code{method} is \code{"num"}, a numerical approximation of
##' the distribution is used.  When \code{method} is \code{"sim"},
##' \code{nSamp} samples of the exponential distribution with the same
##' size as \code{x} are drawn and the LR statistic is computed for
##' each sample. The \eqn{p}-value is simply the estimated probability
##' that a simulated LR is greater than the observed LR. Finally when
##' \code{method} is \code{"asymp"}, the asymptotic distribution is
##' used.
##' 
##' @title Likelihood Ratio test of exponentiality vs. GPD
##'
##' @param x Numeric vector of sample values.
##'
##' @param alternative Character string describing the alternative
##' distribution.
##'
##' @param method Method used to compute the \eqn{p}-value.
##'
##' @param nSamp Number of samples for a simulation, if \code{method}
##' is \code{"sim"}.
##'
##' @return A list of results.
##'
##' @author Yves Deville
##'
##' @note For the Lomax alternative, the distribution of the test
##' statistic has \emph{mixed type}: it can take any positive value as
##' well as the value \eqn{0} with a positive probability mass. The
##' probability mass is the probability that the ML estimate of the
##' GPD shape parameter is negative, and a good approximation of it is
##' provided by the \code{\link{pGreenwood1}} function. Note that this
##' probability converges to its limit \eqn{O.5} \emph{very slowly},
##' which suggests that the asymptotic distribution provides poor
##' results for medium sample sizes, say \eqn{< 100}.
##' 
##' 
LRExp.test <- function(x, 
                       alternative = c("lomax", "GPD", "gpd", "maxlo"),
                       method = c("num", "sim", "asymp"),
                       nSamp = 15000,
                       simW = FALSE){

    eps <- 1e-4
    alternative <- match.arg(alternative)
    method <- match.arg(method)
    
    n <- length(x)
    w <- LRExp(x, alternative = alternative)
    ## cat(sprintf("alternative = %s w = %6.4f\n", alternative, w))
    if ((method == "num") && (n > 500L)) method <- "asymp"
    
    if (method == "num") {
        p <- seq(from = 0.01, to = 0.99, by = 0.01)
        if (alternative == "lomax") {
            pMax <- 1 - pGreenwood1(n)
            if (w < eps) {
                pVal <- pMax
            } else {
                TableL <- qStat(p = p, n = n, type = "logLRLomax")
                nL <- length(TableL$q)
                if (w > TableL$q[nL]) {
                    FL <- 1.0
                } else {
                    FL <- spline(x = TableL$q, y = TableL$p, xout = w)$y
                }
                pVal <- 1 - FL
            }
        } else if (alternative == "maxlo") {
            ## for the maxlo alternative, the sampling distribution is
            ## deduced from that of the GPD case, say FG, and that for
            ## the Lomax case, say FL. Care must be taken to large
            ## values of w since they are very likely to be greater
            ## than the grid max.
            pMax <- pGreenwood1(n)
            if (w < eps) {
                pVal <- pMax
            } else {
                TableG <- qStat(p = p, n = n, type = "logLRGPD")
                nG <- length(TableG$q)
                if (w > TableG$q[nG]) {
                    FG <- 1.0
                } else {
                    FG <- spline(x = TableG$q, y = TableG$p, xout = w)$y
                }
                TableL <- qStat(p = p, n = n, type = "logLRLomax")
                nL <- length(TableL$q)
                if (w > TableL$q[nL]) {
                    FL <- 1.0
                } else {
                    FL <- spline(x = TableL$q, y = TableL$p, xout = w)$y
                }
                pVal <- FL - FG
            }
        } else  if (alternative %in% c("GPD", "gpd")) {
            pMax <- 1
            table <- qStat(p = p, n = n, type = "logLRGPD")
            pVal <- 1 - spline(x = table$q, y = table$p, xout = w)$y
        }
        if (pVal < 0) pVal <- 0
        if (pVal > pMax) pVal <- pMax
    } else if (method == "sim") {
        X <- matrix(rexp(n * nSamp), nrow = nSamp)
        W <- apply(X, 1, LRExp, alternative = alternative)
        pVal <- mean(W > w)
    } else if (method == "asymp") {
        if (alternative %in% c("lomax", "maxlo")) {
            pVal <- 0.5 * pchisq(w, df = 1, lower.tail = FALSE)
        } else  if (alternative %in% c("GPD", "gpd")) {
            pVal <- pchisq(w, df = 1, lower.tail = FALSE)
        }
    }
        
    L <- list(statistic = w, df = n,
              p.value = pVal,
              method = "LR test of exponentiality",
              alternative = alternative)
    
    if (method == "sim" && simW) L[["W"]] <- W
    
    L
    
}

##' Likelihood-Ratio statistic for the Gumbel distribution vs. GEV.
##'
##' The Likelihood-Ratio statistic is actually \eqn{W:=-2 \log
##' \textrm{LR}}{-2 log LR} where LR is the ratio of the likelihoods
##' \emph{Gumbel} to \emph{alternative distribution}. 
##' @title  Likelihood-Ratio statistic
##' @param x Numeric vector of sample values.
##' @param alternative Character string describing the alternative.
##' 
##' @return The LR statistic value.
##' 
##' @author Yves Deville
##'
##' @note When the alternative is \code{"frechet"},
##' the statistic has a distribution of mixed type under the null
##' hypothesis of a Gumbel distribution.
##' 
LRGumbel <- function(x,
                     alternative = c("frechet", "GEV")){

    alternative <- match.arg(alternative)
    n <- length(x)

    fl <- fgev(x, std.err = FALSE)
    if (fl$convergence != "successful") stop("estimation of GEV parameters failed")

    fl0 <- fgev(x, shape = 0, std.err = FALSE)
    if (fl0$convergence != "successful") stop("estimation of Gumbel parameters failed")
   
    if ((alternative == "frechet") && (fl$estimate["shape"] < 0)) {
        W <- 0
    } else {
        W <- fl0$deviance - fl$deviance
    }

    W
}
    
##' Likelihood Ratio test of Gumbel vs. GEV
##'
##' The asymptotic distribution of the Likelihood-ratio statistic is
##' known. For the GEV alternative, this is a chi-square distribution
##' with one df.  For the Frechet alternative, this is the
##' distribution of a product \eqn{XY} where \eqn{X} and \eqn{Y} are
##' two independent random variables following a Bernouilli
##' distribution with probability parameter \eqn{p = 0.5} and a
##' chi-square distribution with one df.
##'
##' When \code{method} is \code{"num"}, a numerical approximation of
##' the distribution is used.  When \code{method} is \code{"sim"},
##' \code{nSamp} samples of the Gumbel distribution with the same size
##' as \code{x} are drawn and the LR statistic is computed for each
##' sample. The \eqn{p}-value is simply the estimated probability that
##' a simulated LR is greater than the observed LR. Finally when
##' \code{method} is \code{"asymp"}, the asymptotic distribution is
##' used.
##' 
##' @title Likelihood Ratio test for the Gumbel distribution
##' 
##' @param x Numeric vector of sample values.
##'
##' @param alternative Character string describing the alternative
##' distribution.
##'
##' @param method Method used to compute the \eqn{p}-value.
##'
##' @param nSamp Number of samples for a simulation, if \code{method}
##' is \code{"sim"}.
##'
##' @return A list of results.
##'
##' @author Yves Deville
##'
##' @note For the Frechet alternative, the distribution of the test statistic has 
##' \emph{mixed type}: it can take any positive value as well as the value \eqn{0}
##' with a positive probability mass. The probability mass is the probability 
##' that the ML estimate of the GEV shape parameter is negative.
##'
##' When \code{method} is \code{"sim"}, the computation can be slow because
##' each of the \code{nSamp} simulated values requires two optimisations.
##' 
##' @examples
##' set.seed(1234)
##' x <- rgumbel(60)
##' res <- LRGumbel.test(x) 
##' 
LRGumbel.test <- function(x, 
                          alternative = c("frechet", "GEV"),
                          method = c("num", "sim", "asymp"),
                          nSamp = 1500,
                          simW = FALSE){
    
    alternative <- match.arg(alternative)
    method <- match.arg(method)

    n <- length(x)
    w <- LRGumbel(x, alternative = alternative)

    if (method == "num") {
        p <- seq(from = 0.01, to = 0.99, by = 0.01)
        if (alternative == "frechet") {
            table <- qStat(p = p, n = n, type = "logLRFrechet")
        } else  if (alternative == "GEV") {
            table <- qStat(p = p, n = n, type = "logLRGEV")
        }
        pVal <- 1 - spline(x = table$q, y = table$p, xout = w)$y
        if (pVal < 0 ) pVal <- 0
        if (pVal > 1 ) pVal <- 1
    } else if (method == "sim") {
        X <- matrix(rgumbel(n * nSamp), nrow = nSamp)
        W <- apply(X, 1, LRGumbel, alternative = alternative)
        pVal <- mean(W > w)
    } else if (method == "asymp") {
        if (alternative == "frechet") {
            pVal <- 0.5 * pchisq(w, df = 1, lower.tail = FALSE)
        } else if (alternative == "GEV") {
            pVal <- pchisq(w, df = 1, lower.tail = FALSE)
        }
    }
    
    L <- list(statistic = w, df = n,
              p.value = pVal,
              method = "LR test for Gumbel",
              alternative = alternative)

    if (method == "sim" && simW) L[["W"]] <- W
    
    L

}
    
