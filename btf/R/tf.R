##' approximate trend filtering via MM algorithm
##'
##' Uses majorization-minimization technique to approximate trend filtering fit.
##' 
##' @param y observed data
##' @param x inputs corresponding to observations 
##' @param k order of fit
##' @param l vector of penalty parameters lambda
##' @param D matrix Delta of order k+1
##' @param eps error adjustment to majorization function
##' @param tau convergence threshold
##' @param max_iter maximum number of iterations allowed
##' @author Edward A. Roualdes
##' @export
tf <- function(y, x=NULL, k=2, l, D=NULL, eps=1e-8, tau=1e-5, max_iter=500) {
    n <- length(y)
    if (missing(l)) stop("Please provide a penalty parameter value l.")
    tf_approx(y, l, D, k, eps=eps, tau=tau, max_iter=max_iter)
}
