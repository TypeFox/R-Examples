#' Kiefer-Wolfowitz estimator for Cosslett (1983) estimator

#' Kiefer-Wolfowitz-Cosslett (1983) estimator for binary response model.  In the
#' primal form of the problem the pseudo log likelihood is:
#' 
#' 	\deqn{l(f|y) =  sum_i [ y_i \log \sum_j (I(v_j <= x_i) * f_j) + 
#' 		(1 - y_i) \log \sum_j (I(v_j > x_i) * f_j) ]}
#'
#' as usual the implementation used here solves the corresponding dual problem.
#' Cumsum of the output y gives the CDF of the unobserved utility difference.
#'
#' @param y is the binary outcome
#' @param x is the observed utility difference between two choices, it would be
#' possible to extend this to make x a linear (index) function of some parameters
#' @param v the unobserved utility difference taking values on a grid 
#' @param weights  replicate weights for x obervations, should sum to 1 
#' @param ... optional parameters to be passed to KWDual to control optimization
#' @return an object of class density with the components:
#' 	\item{x}{points of evaluation of the mixing density} 
#' 	\item{y}{function values of the mixing density at x} 
#'	\item{status}{exit code from the optimizer}
#' @author Jiaying Gu
#' @references Kiefer, J. and J. Wolfowitz Consistency of the Maximum
#' Likelihood Estimator in the Presence of Infinitely Many Incidental
#' Parameters \emph{Ann. Math. Statist}. Volume 27, Number 4 (1956), 887-906.
#' @keywords nonparametric
#' @export

Cosslett <- function(x, y, v = 300, weights = NULL, ...) {
    n <- length(x)
    eps <- 1e-4
    if (length(v) == 1) 
        v <- seq(min(x)-eps, max(x) + eps, length = v)
    m <- length(v)
    if(length(weights)) w <- weights
    else w <- rep(1, n)/n
    d <- diff(v)
    v <- (v[-1] + v[-m])/2
    m <- length(v)
    A <- matrix(0,n,m)
    for (i in 1:n){
    	if (y[i] >0) A[i,] = (x[i]>=v)
    	if (y[i]==0) A[i,] = 1 - (x[i]>=v)
    	}	
    f <- KWDual(A, d, w, ...)
    z <- list(x = v, y = f$f, status = f$status)
    class(z) <- c("Cosslett", "density")
    return(z)
}
