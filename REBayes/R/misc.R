#' Integration by Trapezoidal Rule
#' 
#' Integration by Trapezoidal Rule
#' 
#' Crude Riemann sum approximation.
#' 
#' @param x points of evaluation
#' @param y function values
#' @return A real number.
#' @author R. Koenker
#' @keywords utility
#' @export
traprule <- function(x,y) sum(diff(x) * (y[-1] + y[-length(y)]))/2


#' L1norm for piecewise linear functions
#' 
#' Intended to compute the L1norm of the difference between two distribution
#' functions.
#' 
#' Both F and G should be of class \code{stepfun}, and they should be
#' non-defective distribution functions.  There are some tolerance issues in
#' checking whether both functions are proper distribution functions at the
#' extremes of their support.  For simulations it may be prudent to wrap
#' \code{L1norm} in \code{try}.
#' 
#' @param F A stepfunction
#' @param G Another stepfunction
#' @param eps A tolerance parameter
#' @return A real number.
#' @author R. Koenker
#' @keywords utility
#' @export
#' @importFrom graphics plot.default
#' @importFrom stats approxfun is.stepfun knots rnorm runif var
#' @examples
#' 
#' # Make a random step (distribution) function with Gaussian knots
#' rstep <- function(n){
#'         x <- sort(rnorm(n))
#'         y <- runif(n)
#'         y <- c(0,cumsum(y/sum(y)))
#'         stepfun(x,y)
#'         }
#' F <- rstep(20)
#' G <- rstep(10)
#' S <- L1norm(F,G)
#' plot(F,main = paste("||F - G|| = ", round(S,4)))
#' lines(G,col = 2)
#' 
L1norm <- function(F,G, eps = 1e-6) {
        if(!is.stepfun(F) || !is.stepfun(G)) stop("Both F and G must be stepfun")
        xk <- sort(c(knots(F), knots(G)))
        n <- length(xk)
        if(!all.equal(F(xk[1] - eps), G(xk[1] - eps))) stop("F(x[1]-) != G(x[1]-)")
        if(!all.equal(F(xk[n]), G(xk[n]))) stop("F(x[n]) != G(x[n])")
        dy <- (abs(F(xk) - G(xk)))[-n]
        dx <- diff(xk)
        sum(dy * dx)
        }

