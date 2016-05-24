
##' Cubic Ceschino interpolation.
##'
##' Cubic Ceschino interpolation
##' @title Ceschino cubic interpolation
##' 
##' @param x Numeric vector of design abscissas.
##' 
##' @param y Numeric vector of ordinates.
##'
##' @param xout Numeric vector of new points where interpolation must
##' be done.
##'
##' @param deriv Logical or binary integer. Compute the (first order)
##' derivative?
##' 
##' @return A list with the following elements.
##'
##' \item{x, y}{
##'
##' Coordinates of interpolated points, thus \code{x} is a copy of the
##' input \code{xout}.
##'
##' }
##' \item{CB}{
##'
##' A matrix with \code{length(xout)} rows and \code{length(x)}
##' columns as returned by \code{\link{cardinalBasis_ceschino}}. The
##' columns are the values or \code{xout} of the \eqn{n} Cardinal
##' Basis functions where \eqn{n} is the length of \code{x}.
##' 
##' }
##' \item{deriv}{
##' 
##'  Derivation order used.
##'
##' }
##' \item{knots}{
##'
##'  The knots used in the interpolation. This is a sorted copy of the
##' input \code{x}.
##'
##' }
##' 
##' @author Fortran code by Alain Hebert.
##'
##' @references A. Hebert (2013) \emph{Revisiting the Ceschino
##' Interpolation Method}.
##' \href{http://www.intechopen.com/download/pdf/21941}{link}.
##'
##' 
##' @examples
##' set.seed(12345)
##' n <- 6; nout <- 300L
##' x <- sort(runif(n))
##' xout <- sort(runif(nout, min = x[1], max = x[n]))
##' y <- sin(2 * pi * x)
##' cI0 <- interp_ceschino(x = x, xout = xout, y = y)
##'
##' ## compare with a natural spline
##' require(splines)
##' spI <- interpSpline(x, y)
##' spPred0 <- predict(spI, xout)
##' 
##' plot(xout, sin(2 * pi * xout), type = "l", col = "black", lwd = 2,
##'      xlab = "x", ylab = "f(x)", main = "Interpolations")
##' abline(v = x, col = "gray")
##' lines(cI0, type = "l", col = "SpringGreen3", lty = 2, lwd = 2)
##' lines(spPred0, type = "l", col = "SteelBlue2", lty = 3, lwd = 2)
##' points(x, y, type = "p", pch = 21, col = "red",
##'        bg = "yellow", lwd = 2)
##' legend("topright", legend = c("true", "Ceschino", "nat. spline"),
##'         col = c("black", "SpringGreen3", "SteelBlue2"),
##'         lty = 1:3, lwd = rep(2, 3))
##' 
##' ## derivative estimation
##' cI1 <- interp_ceschino(x =x, xout = xout, y = y, deriv = 1)
##' spPred1 <- predict(spI, xout, deriv = 1)
##' 
##' plot(xout, 2 * pi * cos(2 * pi * xout), type = "l", col = "black", lwd = 2,
##'      xlab = "x", ylab = "fprime(x)", main = "Derivatives")
##' abline(v = x, col = "gray")
##' lines(cI1, type = "l", col = "SpringGreen3", lty = 2, lwd = 2)
##' lines(spPred1, type = "l", col = "SteelBlue2", lty = 3, lwd = 2)
##' points(x, 2 * pi * cos(2 * pi * x), type = "p", pch = 21,
##'        col = "red", bg = "yellow", lwd = 2)
##' legend("bottomright", legend = c("true", "Ceschino", "nat. spline"),
##'        col = c("black", "SpringGreen3", "SteelBlue2"),
##'        lty = 1:3, lwd = rep(2, 3))
##'
interp_ceschino <- function(x, y = NULL, xout, deriv = 0) {

    deriv <- as.integer(deriv)
    if (!(deriv %in% c(0L, 1L))) stop("bad value for 'deriv'")
    
    nout <- length(xout)
    n <- length(x)
    if (length(y) != n) stop("'y' must be a numeric vector with the same length as x")
    
    if (n <= 2L) stop("'n' must be >= 3")
    
    x <- sort(x)
    if ( any( (xout < x[1]) | (xout > x[n]) )) stop("'xout' not in the range of 'x'")
    H <- matrix(0, nrow = nout, ncol = n)
    Work = matrix(0, nrow = 3, ncol = n)
    
    for (i in 1:nout){ 
        
        res <- .Fortran("alterp",
                        1L,                       ## LCUBIC
                        as.integer(n),            ## N
                        as.double(x),             ## X
                        as.double(xout[i]),       ## XNEW
                        as.integer(deriv),        ## DERIV
                        CB = as.double(H[i, ]),   ## Cardinal Basis
                        Work = as.double(Work),   ## Added to avoid allocations
                        PACKAGE = "smint")
        H[i, ] <- res$CB
    }
    
    yout <- H %*% y
    
    list(x = xout, y = yout, CB = H, deriv = deriv, knots = x)
}



                
