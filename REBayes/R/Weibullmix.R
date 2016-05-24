#' NPMLE for Weibull Mixtures
#' 
#' Kiefer-Wolfowitz NPMLE for Weibull Mixtures of scale parameter
#' 
#' Kiefer Wolfowitz NPMLE density estimation for Weibull scale mixtures. The
#' histogram option is intended for relatively large problems, say n > 1000,
#' where reducing the sample size dimension is desirable. By default the grid
#' for the binning is equally spaced on the support of the data.  Parameterization:
#' f(t|alpha, lambda) = alpha * exp(v) * (lambda * t )^(alpha-1) * 
#' exp(-(lambda * t)^alpha * exp(v)); shape = alpha; scale = lambda^(-1) * (exp(v))^(-1/alpha)
#' 
#' @param x Survival times
#' @param v Grid values for mixing distribution
#' @param u Grid values for histogram bins, if needed
#' @param alpha Shape parameter for Weibull distribution
#' @param lambda Scale parameter for Weibull Distribution; must either have length 1, or length
#' equal to \code{length(x)} the latter case accommodates the possibility of a linear predictor
#' @param hist If TRUE aggregate to histogram counts
#' @param weights  replicate weights for x obervations, should sum to 1 
#' @param ...  optional parameters passed to KWDual to control optimization 
#' @return An object of class density with components 
#' 	\item{x}{points of evaluation on the domain of the density} 
#' 	\item{y}{estimated function values at the points x of the mixing density} 
#' 	\item{logLik}{Log likelihood value at the proposed solution} 
#' 	\item{dy}{Bayes Rule estimates of mixing parameter} 
#' 	\item{status}{exit code from the optimizer}
#' @author Roger Koenker and Jiaying Gu
#' @seealso \code{Gompertzmix}
#' @references Kiefer, J. and J. Wolfowitz Consistency of the Maximum
#' Likelihood Estimator in the Presence of Infinitely Many Incidental
#' Parameters \emph{Ann. Math. Statist}. Volume 27, Number 4 (1956), 887-906.
#' @keywords nonparametric
#' @importFrom stats dweibull
#' @export
Weibullmix <- function(x, v = 300, u = 300, alpha, lambda = 1, hist = FALSE, weights = NULL,  ...){
    n <- length(x)
    eps <- 1e-4
    if (length(v) == 1) 
        v <- seq(min(-alpha*log(lambda * x))-eps,max(-alpha*log(lambda * x))+eps,length=v)
	
    w <- weights
    if (hist) {
	if(length(w)) stop("weights not allowed with hist option")
	if(length(u) == 1)
	    u <- seq(min(x) - eps, max(x) + eps, length = u)
	m <- length(u)
        w <- tabulate(findInterval(x, u))
        x <- (u[-1] + u[-m])/2
        wnz <- (w > 0)
        w <- w[wnz]/sum(w[wnz])
        x <- x[wnz]
    }
    if(!length(w)) w <- rep(1,n)/n
    m <- length(v)
    d <- diff(v)
    v <- (v[-1]+v[-m])/2
    if(length(lambda) == 1)
        A <- outer(x, (1/lambda) * (exp(v))^(-1/alpha), FUN = dweibull, shape = alpha)
    else if(length(lambda) == n){
       A <- matrix(0, nrow = n, ncol = length(v))
       for (j in 1:length(v)) {
           s <- (1/lambda) * (exp(v[j]))^(-1/alpha)
           A[, j] <- dweibull(x, shape = alpha, scale = s)
           }
    }
    else
        stop("lambda is of wrong dimension")
    f = KWDual(A, d, w, ...)
    logLik <- n * sum(w * log(f$g))
    dy <- as.vector((A %*% (f$f * d * v))/f$g)  
    z <- list(x = v, y = f$f, g = f$g, logLik = logLik, dy = dy, status = f$status)
    class(z) <- c("Weibullmix", "density")
    return(z)
}
