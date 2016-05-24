
##==========================================================================
##' Jackson's statistic for the exponentiality test.
##'
##' The value(s) of the statistic are the ratio of two weighted means
##' of the order statistics.
##' 
##' @title Jackson's statistic
##'
##' @param x Numeric vector or matrix. In the second case, rows are
##' considered as samples.
##'
##' @param norm Logical: if \code{TRUE}, the statistic is normalized
##' by using the \emph{asymptotic} mean and standard deviation,
##' respectively 2 and 1.
##'
##' @return A numeric vector of length \code{1} when \code{x} is a vector,
##' or with length \code{nrow(x)} when \code{x} is a matrix. 
##'
##' @references 
##' J. Beirlant and T. de Weit and Y. Goegebeur(2006)
##' A Goodness-of-fit Statistic for Pareto-Type Behaviour,
##' \emph{J. Comp. Appl. Math.}, 186(1), pp. 99-116}
##'
Jackson <- function(x, norm = FALSE) {

    if (is.null(dim(x))) {
       x <- sort(x)
       t <- cumsum(1 / rev(seq(along = x)))
       J <- sum(t * x) / sum(x)
       if (norm) J <- sqrt(length(x)) * (J -2)
       return(J)
    } else {
        if (!is.matrix(x)) {
            stop("'x' must be a vector or a matrix")
        }
        n <- ncol(x)
        s <- apply(x, 1, sum)
        x <- t(apply(x, 1, sort))
        t <- cumsum(1 / (n:1L))
        x <- sweep(x, MARGIN = 2, STATS = t, FUN = "*")
        tx <- apply(x, 1, sum)
        J <- tx / s
        if (norm) J <- sqrt(n) * (J -2)
        return(J)
    }
    
}

##==========================================================================
##' Random generation of Jackson's statistic.
##'
##' The values are simulated using the expression of the
##' statistic involving the normalised spacings, so that no sorting
##' operation is required. Hence the computation is faster than generating
##' exponential samples and then computing their CV2.
##' @title Random generation of Jacskon's statistic
##' @param n number of simulated values.
##' @param size size of the sample.
##' @return a numeric vector of length \code{n}.
##' @author Yves Deville
##' @examples
##' set.seed(234)
##' J <- rJackson(n = 100, size = 30)
##' mean(J)  ## nearly 2
##' sqrt(30) * sd(J)    ## nearly 1
rJackson <- function(n , size) {
    w <- c(0, Renext::Hpoints(size))
    w <- 1 + w[1L:size]
    ## use normalised spacings
    NS <- matrix(rexp(n * size), nrow = n, ncol = size)
    Den <- apply(NS, 1, sum)
    Num <- sweep(NS, MARGIN = 2, STATS = w, FUN = "*")
    Num <- apply(Num, 1, sum)
    Num / Den
}

##==========================================================================
##' Jackson's test of exponentiality
##'
##' Compute the Jackson's test of exponentiality.  The test statistic
##' is the ratio of weighted sums of the order statistics. Both sums
##' can also be written as weighted sums of the scalings.
##' 
##' The Jackson's statistic for a sample of size \eqn{n} of the exponential 
##' distribution can be shown to be approximatively normal. More precisely 
##' \eqn{\sqrt{n}(J_n -2)}{sqrt(n)*(Jn -2)} has approximatively a standard normal
##' distribution. This distribution is used to compute the \eqn{p}-value when 
##' \code{method} is \code{"asymp"}. When \code{method} is \code{"num"}, a numerical
##' approximation of the distribution is used. Finally, when \code{method} is 
##' \code{"sim"} the \eqn{p}-value is computed by simulating \code{nSamp} samples
##' of size \code{length(x)} and estimating the probability to have a Jackson's statistic
##' larger than that of the 'observed' \code{x}. 
##' @title Jackson's test of exponentiality
##' @param x numeric vector or matrix.
##' @param method Character: choice of the method used to compute the
##' \eqn{p}-value. See the \bold{Details} section.
##' @param nSamp Number of samples used to compute the \eqn{p}-value
##' if \code{method} is \code{"sim"}. 
##' @return A list of results
##' @author Yves Deville
##' @seealso The \code{\link{LRExp.test}} function.
##' @note Jackson's test of exponentiality works fine for a Lomax alternative
##' (GPD whith heavy tail). It then reaches nearly the same power as a Likelihood Ratio
##' test, while being easier to implement.
Jackson.test <- function(x,
                         method = c("num", "sim", "asymp"),
                         nSamp = 15000){
    
    method <- match.arg(method)
    n <- length(x)
    j <- Jackson(x)
    
    if (method == "num") {
        p <- seq(from = 0.01, to = 0.99, by = 0.01)
        table <- qStat(p = p, n = n, type = "Jackson")
        pVal <- 1 - spline(x = table$q, y = table$p, xout = j)$y
        if (pVal < 0 ) pVal <- 0
        if (pVal > 1 ) pVal <- 1
    } else if (method == "sim") {
        J <- rJackson(nSamp, size = n)
        pVal <- mean(J > j)
    } else if (method == "asymp") {
        jMod <- sqrt(n) * (j - 2) 
        pVal <- pnorm(jMod, lower.tail = FALSE)
    }
    pVal <- round(pVal, digits = 3)
    
    list(statistic = j, df = n,
         p.value = pVal,
         method = "Jackson's test of exponentiality")
        
}
