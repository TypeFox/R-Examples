#' NPMLE for Gompertz Mixtures
#' 
#' Kiefer-Wolfowitz NPMLE for Gompertz Mixtures of scale parameter
#' 
#' Kiefer Wolfowitz NPMLE density estimation for Gompertz scale mixtures. The
#' histogram option is intended for relatively large problems, say n > 1000,
#' where reducing the sample size dimension is desirable. By default the grid
#' for the binning is equally spaced on the support of the data. 
#' Parameterization: f(t|alpha,theta,v) = theta * exp(v) * exp(alpha * t) * 
#'          exp(-(theta/alpha) * exp(v) * (exp(alpha*t)-1))
#' 
#' @param x Survival times
#' @param v Grid values for mixing distribution
#' @param u Grid values for mixing distribution
#' @param alpha Shape parameter for Gompertz distribution
#' @param theta Scale parameter for Gompertz Distribution
#' @param hist If TRUE aggregate to histogram counts
#' @param weights  replicate weights for x obervations, should sum to 1 
#' @param ... optional parameters passed to KWDual to control optimization 
#' @return An object of class density with components 
#' 	\item{x}{points of evaluation on the domain of the density} 
#' 	\item{y}{estimated function values at the points x, the mixing density} 
#' 	\item{logLik}{Log likelihood value at the proposed solution} 
#' 	\item{dy}{Bayes rule estimates of theta at observed x} 
#' 	\item{status}{exit code from the optimizer}
#' @author Roger Koenker and Jiaying Gu
#' @seealso \code{Weibullmix}
#' @references Kiefer, J. and J. Wolfowitz Consistency of the Maximum
#' Likelihood Estimator in the Presence of Infinitely Many Incidental
#' Parameters \emph{Ann. Math. Statist}. Volume 27, Number 4 (1956), 887-906.
#' @keywords nonparametric
#' @export
Gompertzmix <- function(x, v = 300, u = 300, alpha, theta, hist = FALSE, weights = NULL, ...){
    requireNamespace("reliaR", quietly = TRUE)
    n <- length(x)
    eps <- 1e-4
    if (length(v) == 1) # Support is bounded by max and min of MLE for v
        v <- seq(min(-log((theta/alpha)*(exp(alpha*x)-1)))-eps,
	   max(-log((theta/alpha)*(exp(alpha*x)-1)))+eps,length=v) 
    m <- length(v)
    if (hist) {
        if(length(u) == 1) 
	    u <- seq(min(x) - eps, max(x) + eps, length = u)
	m <- length(u)
        w <- tabulate(findInterval(x, u))
        x <- (u[-1] + u[-m])/2
        wnz <- (w > 0)
        w <- w[wnz]/sum(w[wnz])
        x <- x[wnz]
    }
    if(length(weights)) w <- weights
    else w <- rep(1,n)/n
    d <- diff(v)
    v <- (v[-1] + v[-m])/2
    A <- outer(x,theta *exp(v),FUN=reliaR::dgompertz,alpha=alpha)
    f = KWDual(A, d, w, ...)
    logLik <- n * sum(w * log(f$g))
    dy <- as.vector((A %*% (f$f * d * v))/f$g)  # Bayes rule for v
    z <- list(x = v, y = f$f, g = f$g, logLik = logLik, dy = dy, status = f$status)
    class(z) <- c("Gompertzmix", "density")
    return(z)
}
