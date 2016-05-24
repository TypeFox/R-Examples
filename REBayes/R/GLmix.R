#' Kiefer-Wolfowitz NPMLE for Gaussian Location Mixtures
#' 
#' Kiefer Wolfowitz Nonparametric MLE for Gaussian Location Mixtures
#' 
#' Kiefer Wolfowitz MLE as proposed by Jiang and Zhang for
#' the Gaussian compound decision problem.  The histogram option is intended
#' for large problems, say n > 1000, where reducing the sample size dimension
#' is desirable. When \code{sigma} is heterogeneous and \code{hist = TRUE} the
#' procedure tries to do separate histogram binning for distinct values of
#' \code{sigma}, however this is only feasible when there are only a small
#' number of distinct \code{sigma}. By default the grid for the binning is
#' equally spaced on the support of the data. This function does the normal
#' convolution problem, for gamma mixtures of variances see \code{GVmix}, or
#' for mixtures of both means and variances \code{TLVmix}.  
#' 
#' The predict method for \code{GLmix} objects will compute means, medians or
#' modes of the posterior according to whether the \code{Loss} argument is 2, 1
#' or 0.
#' 
#' @param x Data: Sample Observations
#' @param v Undata: Grid Values defaults equal spacing of with v bins, when v is
#' a scalar
#' @param sigma scale parameter of the Gaussian noise, may take vector values
#' of length(x)
#' @param hist If TRUE then aggregate x to histogram bins, when sigma is vector
#' valued this option is inappropriate unless there are only a small number of
#' distinct sigma values.
#' @param histm histogram bin boundaries, equally spacing with \code{histm} 
#' bins when  scalar.
#' @param weights  replicate weights for x obervations, should sum to 1 
#' @param ... other parameters to pass to KWDual to control optimization
#' @return An object of class density with components: 
#' 	\item{x}{points of  evaluation on the domain of the density} 
#' 	\item{y}{estimated function values at the points v, the mixing density} 
#' 	\item{g}{the estimated mixture density function values at x} 
#' 	\item{logLik}{Log likelihood value at the proposed solution} 
#' 	\item{dy}{prediction of mean parameters for each observed x value via Bayes Rule} 
#' 	\item{status}{exit code from the optimizer}
#' @author Roger Koenker
#' @references Kiefer, J. and J. Wolfowitz Consistency of the Maximum
#' Likelihood Estimator in the Presence of Infinitely Many Incidental
#' Parameters \emph{Ann. Math. Statist}.  Volume 27, Number 4 (1956), 887-906.
#' 
#' Jiang, Wenhua and Cun-Hui Zhang General maximum likelihood empirical Bayes
#' estimation of normal means \emph{Ann. Statist.}, Volume 37, Number 4 (2009),
#' 1647-1684.
#'
#' Koenker, R and I. Mizera, (2013) ``Convex Optimization, Shape Constraints,
#' Compound Decisions, and Empirical Bayes Rules,'' \emph{JASA}, 109, 674--685.
#' @keywords nonparametric
#' @importFrom stats dnorm
#' @export
GLmix <- function (x, v = 300, sigma = 1, hist = FALSE, histm = 300, weights = NULL, ...)
{
    n <- length(x)
    eps <- 1e-4
    if (length(v) ==1) 
        v <- seq(min(x) - eps, max(x) + eps, length = v)
    m <- length(v)
    weights <- weights/sum(weights)
    w <- weights
    if (hist) {
      histbin <- function(x, m = histm, eps = 1e-06, weights = weights) {
        u <- seq(min(x) - eps, max(x) + eps, length = m)
        xu <- findInterval(x, u)
	txu <- tabulate(xu)
        midu <- (u[-1] + u[-m])/2
        wnz <- (txu > 0)
	if(length(weights)){ 
	    if(length(weights) == length(x)) 
		w <- as.numeric(tapply(weights, xu, sum))
	    else
		stop("length(weights) not equal length(x)")
	    }
	else
	    w <- txu[wnz]/sum(txu[wnz])
        list(midu = midu[wnz], w = w)
        }
      if(length(sigma) == 1){
	  h <- histbin(x, m = histm, weights = weights)
	  x <- h$midu
	  w <- h$w
	}
      else { # create sigma bins for histogram
	    sus <- sort(unique(sigma))
	    us <- match(sigma, sus)
	    nus <- table(us)
	    if(min(nus) < 100) stop("too few obs in some sigma bin")
	    h <- as.list(1:length(nus))
	    for(i in 1:length(sus)){
		if(length(weights))
		    h[[i]] <- histbin(x[us == i],m = histm, weights = weights[us == i])
		else
		    h[[i]] <- histbin(x[us == i],m = histm, weights = weights)
	    }
	    x <- unlist(lapply(h, function(f) f$midu))
	    w <- unlist(lapply(h, function(f) f$w))
	    w <- w/sum(w)
	    sigma <- rep(sus,unlist(lapply(h, function(f) length(f$midu))))
      }
    }
    if(!length(w)) w <- rep(1, n)/n
    d <- diff(v)
    v <- (v[-1] + v[-m])/2
    A <- dnorm(outer(x, v, "-"), sd = sigma)
    f <- KWDual(A, d, w, ...)
    y <- f$f
    g <- f$g
    logLik <- n * sum(w * log(g))
    dy <- as.vector((A %*% (y * d * v))/g)
    z <- list(x = v, y = y, g = g, logLik = logLik, 
	sigma = sigma, dx = x, dy = dy, status = f$status)
    class(z) <- c("GLmix", "density")
    return(z)
}
