#' NPMLE for Student t location mixtures
#' 
#' Kiefer Wolfowitz NPMLE for Student t location mixtures
#' 
#' Kiefer Wolfowitz MLE density estimation as proposed by Jiang and Zhang for
#' a Student t compound decision problem.  The histogram option is intended
#' for large problems, say n > 1000, where reducing the sample size dimension
#' is desirable. By default the grid for the binning is equally spaced on the
#' support of the data. Equal spaced binning is problematic for Cauchy data.
#' 
#' @param x Data: Sample Observations
#' @param v bin boundaries defaults to equal spacing of length v
#' @param u bin boundaries for histogram binning: defaults to equal spacing 
#' @param df Number of degrees of freedom of Student base density
#' @param hist If TRUE then aggregate x to histogram weights
#' @param weights  replicate weights for x obervations, should sum to 1 
#' @param ... optional parameters passed to KWDual to control optimization
#' @return An object of class density with components: 
#' 	\item{x}{midpoints of evaluation on the domain of the mixing density} 
#' 	\item{y}{estimated function values at the points x of the mixing density} 
#' 	\item{logLik}{Log likelihood value at the proposed solution} 
#' 	\item{dy}{Bayes rule estimates of location at x} 
#' 	\item{status}{Mosek exit code}
#' @author Roger Koenker
#' @references Kiefer, J. and J. Wolfowitz Consistency of the Maximum
#' Likelihood Estimator in the Presence of Infinitely Many Incidental
#' Parameters \emph{Ann. Math. Statist}.  27, (1956), 887-906.
#' 
#' Jiang, Wenhua and Cun-Hui Zhang General maximum likelihood empirical Bayes
#' estimation of normal means \emph{Ann. Statist.}, 37, (2009), 1647-1684.
#' @seealso GLmix for Gaussian version
#' @keywords nonparametric
#' @importFrom stats dt
#' @export
TLmix <- function(x, v = 300, u = 300, df = 1, hist = FALSE, weights = NULL, ...){

   n <- length(x)
   eps <- 1e-4
   if(length(v) == 1) v <- seq(min(x)-eps, max(x)+eps, length = v)
   if(hist){
      if(length(u) == 1) u <- seq(min(x)-eps,max(x)+eps,length = u)
      w <- tabulate(findInterval(x,u))
      x <- (u[-1] + u[-m])/2
      wnz <- (w > 0)
      w <- w[wnz]/sum(w[wnz])
      x <- x[wnz]
      }
   if(length(weights)) w <- weights
   else w <- rep(1,n)/n
   m <- length(v)
   d <- diff(v)
   v <- (v[-1] + v[-m])/2
   A <- dt(outer(x,v,"-"),df = df) 
   f = KWDual(A, d, w, ...)
   logLik <- n * sum(w * log(f$g))
   dy <- as.vector((A %*% (f$f * d * v))/f$g)
   z <- list(x = v, y = f$f, g = f$g, logLik = logLik, dy = dy, status = f$status)
   class(z) <- c("TLmix", "density")
return(z)
}
