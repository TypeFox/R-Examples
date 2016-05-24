#' NPMLE for Student t non-centrality parameter mixtures
#' 
#' Kiefer Wolfowitz NPMLE for Student t non-centrality parameter mixtures
#' Model: \eqn{y_{ig} = mu_{g} + e_{ig}, e_{ig} ~ N(0,sigma_{g}^{2})}
#' x is the vector of t statistics for all groups, which follows t dist
#' if \eqn{mu_g = 0}, and noncentral t dist if \eqn{mu_g \neq 0}, 
#' with \eqn{ncp_{g} = \mu_g / \sigma_{g}}.
#' This leads to a mixture of t distribution with ncp as the mixing parameter.  
#' df (degree of freedom) is determined by the group size in the simplest case.
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
#' 	\item{g}{estimated function values at the observed points of mixture density} 
#' 	\item{logLik}{Log likelihood value at the proposed solution} 
#' 	\item{dy}{Bayes rule estimates of location at x} 
#' 	\item{status}{Mosek exit code}
#' @author Roger Koenker
#' @references Kiefer, J. and J. Wolfowitz Consistency of the Maximum
#' Likelihood Estimator in the Presence of Infinitely Many Incidental
#' Parameters \emph{Ann. Math. Statist}.  27, (1956), 887-906.
#' 
#' @seealso GLmix for Gaussian version
#' @keywords nonparametric
#' @export


Tncpmix <- function (x, v = 300, u = 300, df = 1,  hist = FALSE, weights = NULL, ...) {
    n <- length(x)
    eps <- 1e-4
    if (length(df)==1) df = rep(df,n)
    if (missing(v)) 
        v <- seq(min(x) - eps, max(x) + eps, length = v)
    if (hist) {
        if (length(u) == 1) 
            u <- seq(min(x) - eps, max(x) + eps, length = u)
	m <- length(u)
        w <- tabulate(findInterval(x, u))
        x <- (u[-1] + u[-m])/2
        wnz <- (w > 0)
        w <- w[wnz]/sum(w[wnz])
        x <- x[wnz]
    }
    if(length(weights)) w <- weights
    else w <- rep(1, n)/n
    m <- length(v)
    d <- diff(v)
    v <- (v[-1] + v[-m])/2
    A = outer(x, v, FUN = dt, df = df)
    f = KWDual(A, d, w, ...)
    y <- f$f
    g <- f$g
    logLik <- n * sum(w * log(g))
    dy <- as.vector((A %*% ( y * d * v))/g) # Bayes rule on noncentrality parameter
    z <- list(x = v, y = y, g = f$g, logLik = logLik, dy = dy,
        status = f$status)
    class(z) <- "density"
    return(z)
}

