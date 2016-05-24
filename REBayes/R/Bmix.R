#' Binomial mixture estimation via Kiefer Wolfowitz MLE
#' 
#' Interior point solution of Kiefer-Wolfowitz NPMLE for mixture of binomials
#' 
#' @param x Count of "successes" for binomial observations
#' @param k Number of trials for binomial observations
#' @param v Grid Values for the mixing distribution defaults to equal
#' spacing of length v on [eps, 1- eps], if v is scalar.
#' @param collapse Collapse observations into cell counts.
#' @param weights  replicate weights for x obervations, should sum to 1 
#' @param ... Other arguments to be passed to KWDual to control optimization
#' @return An object of class density with components: 
#' 	\item{x}{grid midpoints of evaluation of the mixing density} 
#' 	\item{y}{function values of the mixing density at x} 
#' 	\item{g}{estimates of the mixture density at the distinct data values} 
#' 	\item{logLik}{Log Likelihood value at the estimate}
#' 	\item{dy}{Bayes rule estimates of binomial probabilities for distinct data values}
#' 	\item{status}{exit code from the optimizer}
#' @author R. Koenker
#' @references Kiefer, J. and J. Wolfowitz Consistency of the Maximum
#' Likelihood Estimator in the Presence of Infinitely Many Incidental
#' Parameters \emph{Ann. Math. Statist}. 27, (1956), 887-906.
#'
#' Koenker, R and I. Mizera, (2013) ``Convex Optimization, Shape Constraints,
#' Compound Decisions, and Empirical Bayes Rules,'' \emph{JASA}, 109, 674--685.
#' @keywords nonparametric
#' @importFrom stats dbinom
#' @export
Bmix <- function(x, k, v = 300, collapse = TRUE, weights = NULL, ...){

    n <- length(x)
    w <- weights
    if(collapse){ #collapse observations into cell counts
      T <- table(x,k)
      x <- rep(as.numeric(dimnames(T)[[1]]), NCOL(T))
      k <- rep(as.numeric(dimnames(T)[[2]]), each = NROW(T))
      y <- c(T)
      s <- y > 0
      y <- y[s]
      x <- x[s]
      k <- k[s]
      w <- y/sum(y)
      }
    if(!length(w)) w <-  rep(1,n)/n
   eps <- 1e-4
   if(length(v) == 1) v <- seq(eps, 1 - eps, length = v)
   m <- length(v)
   d <- diff(v)
   v <- (v[-1] + v[-m])/2
   A <- outer(x,v,function(x, v, k) dbinom(x,size = k, prob = v), k = k)
   z <- KWDual(A, d, w, ...)
   g <- z$g
   logLik <- n * sum(w * log(g))
   dy <- as.vector((A%*%(z$f * d * v))/g)
   z <- list(x = v, y = z$f, g = g, logLik = logLik, dy = dy, status= z$status)
class(z) <- c("Bmix", "density")
return(z)
}
