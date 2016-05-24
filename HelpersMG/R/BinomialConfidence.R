#' @title Confidence Intervals for Binomial Probabilities
#' @author Rollin Brant, Modified by Frank Harrell and Brad Biggerstaff 
#' @return a matrix or data.frame containing the computed intervals and, optionally, x and n.
#' @param x Vector containing the number of "successes" for binomial variates
#' @param n Vector containing the numbers of corresponding observations
#' @param alpha Probability of a type I error, so confidence coefficient = 1-alpha
#' @param method Character string specifing which method to use. The "all" method only works when x and n are length 1. The "exact" method uses the F distribution to compute exact (based on the binomial cdf) intervals; the "wilson" interval is score-test-based; and the "asymptotic" is the text-book, asymptotic normal interval. Following Agresti and Coull, the Wilson interval is to be preferred and so is the default.
#' @param include.x Logical flag to indicate whether x should be included in the returned matrix or data frame
#' @param include.n Logical flag to indicate whether n should be included in the returned matrix or data frame
#' @param return.df Logical flag to indicate that a data frame rather than a matrix be returned
#' @description Produces 1-alpha confidence intervals for binomial probabilities.\cr
#' @examples
#' HelpersMG:::.BinomialConfidence(0:10,10,include.x=TRUE,include.n=TRUE)
#' HelpersMG:::.BinomialConfidence(46,50,method="all")
#' @export

.BinomialConfidence <- 
function (x, n, alpha = 0.05, method = c("wilson", "exact", "asymptotic", 
    "all"), include.x = FALSE, include.n = FALSE, return.df = FALSE) 
{
    method <- match.arg(method)
    bc <- function(x, n, alpha, method) {
        nu1 <- 2 * (n - x + 1)
        nu2 <- 2 * x
        ll <- if (x > 0) 
            x/(x + qf(1 - alpha/2, nu1, nu2) * (n - x + 1))
        else 0
        nu1p <- nu2 + 2
        nu2p <- nu1 - 2
        pp <- if (x < n) 
            qf(1 - alpha/2, nu1p, nu2p)
        else 1
        ul <- ((x + 1) * pp)/(n - x + (x + 1) * pp)
        zcrit <- -qnorm(alpha/2)
        z2 <- zcrit * zcrit
        p <- x/n
        cl <- (p + z2/2/n + c(-1, 1) * zcrit * sqrt((p * (1 - 
            p) + z2/4/n)/n))/(1 + z2/n)
        if (x == 1) 
            cl[1] <- -log(1 - alpha)/n
        if (x == (n - 1)) 
            cl[2] <- 1 + log(1 - alpha)/n
        asymp.lcl <- x/n - qnorm(1 - alpha/2) * sqrt(((x/n) * 
            (1 - x/n))/n)
        asymp.ucl <- x/n + qnorm(1 - alpha/2) * sqrt(((x/n) * 
            (1 - x/n))/n)
        res <- rbind(c(ll, ul), cl, c(asymp.lcl, asymp.ucl))
        res <- cbind(rep(x/n, 3), res)
        switch(method, wilson = res[2, ], exact = res[1, ], asymptotic = res[3, 
            ], all = res, res)
    }
    if ((length(x) != length(n)) & length(x) == 1) 
        x <- rep(x, length(n))
    if ((length(x) != length(n)) & length(n) == 1) 
        n <- rep(n, length(x))
    if ((length(x) > 1 | length(n) > 1) & method == "all") {
        method <- "wilson"
        warning("method=all will not work with vectors...setting method to wilson")
    }
    if (method == "all" & length(x) == 1 & length(n) == 1) {
        mat <- bc(x, n, alpha, method)
        dimnames(mat) <- list(c("Exact", "Wilson", "Asymptotic"), 
            c("PointEst", "Lower", "Upper"))
        if (include.n) 
            mat <- cbind(N = n, mat)
        if (include.x) 
            mat <- cbind(X = x, mat)
        if (return.df) 
            mat <- as.data.frame(mat)
        return(mat)
    }
    mat <- matrix(ncol = 3, nrow = length(x))
    for (i in 1:length(x)) mat[i, ] <- bc(x[i], n[i], alpha = alpha, 
        method = method)
    dimnames(mat) <- list(rep("", dim(mat)[1]), c("PointEst", 
        "Lower", "Upper"))
    if (include.n) 
        mat <- cbind(N = n, mat)
    if (include.x) 
        mat <- cbind(X = x, mat)
    if (return.df) 
        mat <- as.data.frame(mat, row.names = NULL)
    mat
}
