##' plot btf object
##'
##' Plots a btf object with optional credible intervals.
##' @param x btf object
##' @param t domain of function
##' @param burn size of burn-in,
##' @param est function specifying how draws from the posterior are summarized
##' @param probs numeric 2-vector of credible interval probabilities; if FALSE, no credible intervals are plot
##' @param ... extra arguments
##' @author Edward A. Roualdes
##' @aliases plot.btf
##' @export
plot.btf <- function(x, t=NULL, burn = 1e3, est=median,
                     probs = c(0.025, 0.975), ...) {
    btf <- x
    y <- attr(btf, 'y')
    n <- length(y)
        if (missing(t)) {
        t <- (1:n)/n
    }

    est <- match.fun(est)
    beta <- getPostEst(btf, 'beta', est=est) # point estimates
    q <- getPostEst(btf, 'beta', est=function(z) quantile(z, probs)) # ci
    ry <- range(y); rq <- range(q)
    
    plot(y~t, xlim = range(t), ylim=c(min(ry[1], rq[1]), max(ry[2], rq[2])),
         xlab = 'x', ylab='f(x)', ...)
    lines(beta~t, col=cols['blue2'], lwd=2)
    if ( all(as.logical(probs)) ) {
        lines(q[1,]~t, col = cols['blue2'], lty=2, lwd=2)
        lines(q[2,]~t, col = cols['blue2'], lty=2, lwd=2)        
    }
}
