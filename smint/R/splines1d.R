
##' Cardinal Basis for cubic Ceschino interpolation.
##'
##' This is a simple and raw interface to \code{alterp} Fortran
##' subroutine.
##' 
##' @title Cardinal Basis for cubic Ceschino interpolation
##'
##' @param x Numeric vector of design points.
##'
##' @param xout Numeric vector giving new points.
##'
##' @param cubic Logical. Use cubic interpolation or basic linear?
##'
##' @param deriv Integer or logical. Compute the derivative?
##' 
##' @return A list with the following elements
##'
##' \item{x}{
##'
##' Numeric vector of abscissas at which the basis is evaluated. This
##' is a copy of \code{xout}.
##'
##' }
##' \item{CB}{
##'
##' Matrix of the Cardinal Basis function values.
##'
##' }
##' \item{deriv, cubic}{
##'
##' Copy of input.
##'
##' }
##' \item{method}{
##'
##' Character description of the method involved in the CB determination.
##'
##' }
##'
##' @author Alain Hebert for Fortran code.
##'
##' Yves Deville for R interface.
##'
##' @note This function does not allow extrapolation, so an error will
##' result when \code{xout} contains element outside of the range of
##' \code{x}.
##'
##' @seealso \code{\link{interp_ceschino}} for the related interpolation
##' function, \code{\link{cardinalBasis_natSpline}} and \code{\link{cardinalBasis_lagrange}}
##' for other Cardinal Basis constructions.
##' 
##' @examples
##' set.seed(123)
##' n <- 16L; nout <- 300L
##' x <- sort(runif(n))
##'
##' ## let 'xout' contain n + nout points including nodes 
##' xout <- sort(c(x, runif(nout, min = x[1], max = x[n])))
##' y <- sin(2 * pi * x)
##' res  <- cardinalBasis_ceschino(x, xout = xout, deriv = 0)
##' 
##' matplot(res$x, res$CB, type = "n", main = "Cardinal Basis")
##' abline(v = x, h = 1.0, col = "gray")
##' points(x = x, y = rep(0, n), pch = 21, col = "black",
##'        lwd = 2, bg = "white")
##' matlines(res$x, res$CB, type = "l")
##' 
##' ## interpolation error should be fairly small
##' max(abs(sin(2 * pi * xout) - res$CB \%*\% y))
##' 
cardinalBasis_ceschino <- function(x, xout, cubic = TRUE, deriv = 0) {
    
    nout <- length(xout)
    n <- length(x)
    
    if (!cubic && (n <= 1L)) stop("'n' must be >= 2")
    if (cubic && (n <= 2L)) stop("'n' must be >= 3 when 'cubic' is TRUE")
    
    x <- sort(x)
    if ( any( (xout < x[1]) | (xout > x[n]) )) {
        stop("'xout' not in the range of 'x'")
    }
    H <- matrix(0, nrow = nout, ncol = n)
    Work <- matrix(0, nrow = 3, ncol = n)
    
    for (i in 1:nout){ 
        
        res <- .Fortran("alterp",
                        as.integer(cubic),        ## CUBIC
                        as.integer(n),            ## N
                        as.double(x),             ## X
                        as.double(xout[i]),       ## XNEW
                        as.integer(deriv),        ## DERIV
                        CB = as.double(H[i, ]),   ## Cardinal Basis
                        Work = as.double(Work),   ## Added to avoid allocations
                        PACKAGE = "smint")
        H[i, ] <- res$CB
    }
    
    list(x = xout,
         CB = H,
         deriv = deriv,
         cubic = cubic,
         method = "Ceschino interpolation")
  
}

##' Cardinal Basis for natural cubic spline interpolation.
##'
##' This is a simple and raw interface to \code{splinterp} Fortran
##' subroutine.
##'
##' @title Cardinal Basis for natural cubic spline interpolation
##'
##' @param x Numeric vector of design points.
##'
##' @param xout Numeric vector of new points.
##'
##' @param deriv Integer. Order of derivation. Can be \code{0}, \code{1}
##' or \code{2}.
##' 
##' @return A list with several elements
##'
##' \item{x}{
##'
##' Numeric vector of abscissas at which the basis is evaluated. This
##' is a copy of \code{xout}.
##'
##' }
##' \item{CB}{
##'
##' Matrix of the Cardinal Basis function values.
##'
##' }
##' \item{deriv}{
##'
##' Order of derivation as given on input.
##'
##' }
##' \item{method}{
##'
##' Character description of the method involved in the CB determination.
##'
##' }
##' @author Yves Deville
##'
##' @examples
##' set.seed(123)
##' n <- 16; nout <- 360
##' x <- sort(runif(n))
##' 
##' ##' ## let 'xout' contain n + nout points including nodes 
##' xout <- sort(c(x, seq(from =  x[1] + 1e-8, to = x[n] - 1e-8, length.out = nout)))
##' res  <- cardinalBasis_natSpline(x, xout = xout)
##' 
##' matplot(res$x, res$CB, type = "n", main = "Cardinal Basis")
##' abline(v = x, h = 1.0, col = "gray")
##' points(x = x, y = rep(0, n), pch = 21, col = "black", lwd = 2, bg = "white")
##' matlines(res$x, res$CB, type = "l")
##'
##' ## compare with 'splines'
##' require(splines)
##' y <- sin(2* pi * x)
##' sp <- interpSpline(x, y)
##' test <- rep(NA, 3)
##' der <- 0:2
##' names(test) <- nms <- paste("deriv. ", der, sep = "")
##' for (i in seq(along = der)) {
##'    resDer <- cardinalBasis_natSpline(x, xout = xout, deriv = der[i])
##'    test[nms[i]] = max(abs(predict(sp, xout, deriv = der[i])$y - resDer$CB \%*\% y))
##' }
##' test
##' ## Lebesgue's function
##' plot(x = xout, y = apply(res$CB, 1, function(x) sum(abs(x))), type = "l",
##'      lwd = 2, col = "orangered", main = "Lebesgue\'s function", log = "y",
##'      xlab = "x", ylab = "L(x)")
##' points(x = x, y = rep(1, n), pch = 21, col = "black", lwd = 2, bg = "white")

cardinalBasis_natSpline <- function(x, xout, deriv = 0) {
    
    nout <- length(xout)
    n <- length(x)
    if (n <= 2L) stop("'n' must be >= 3")
    x <- sort(x)
    
    if ( any( (xout < x[1]) | (xout > x[n]) )) stop("'xout' not in the range of 'x'")
    H <- matrix(0, nrow = nout, ncol = n)
    Work <- matrix(0, nrow = 2, ncol = n - 2)
    B <- matrix(0, nrow = n, ncol = n)
    
    for (i in 1L:nout){ 
        
        res <- .Fortran("splinterp",
                        1L,                        ## CUBIC
                        as.integer(n),             ## N
                        as.double(x),              ## X
                        as.double(xout[i]),        ## XNEW
                        as.integer(deriv),         ## DERIV
                        CB = as.double(H[i, ]),    ## Cardinal Basis
                        Work = as.double(Work),    ## to avoid allocations
                        B = as.double(B),          ## idem
                        PACKAGE = "smint")
        H[i, ] <- res$CB

    }
    
    list(x = xout,
         CB = H,
         deriv = deriv,
         knots = xout,
         method = "natural spline")
    
}

##' Cardinal Basis for Lagrange interpolation.
##'
##' This is a simple and raw interface to \code{alterp} Fortran
##' subroutine. It is a wrapper for \code{\link{cardinalBasis_ceschino}}
##' function with \code{cubic = FALSE} and \code{deriv = 0L}.
##' 
##' @title Cardinal Basis for Lagrange (broken line) interpolation
##'
##' @param x Numeric vector of design points.
##'
##' @param xout Numeric vector giving new points.
##' 
##' @return A list with the following elements
##'
##' \item{x}{
##'
##' Numeric vector of abscissas at which the basis is evaluated. This
##' is a copy of \code{xout}.
##'
##' }
##' \item{CB}{
##'
##' Matrix of the Cardinal Basis function values.
##'
##' }
##'
##' @note This function does not allow extrapolation, so an error will
##' result when \code{xout} contains element outside of the range of
##' \code{x}. The function used here is a spline of degree 1 (order 2).
##'
##' @examples
##' set.seed(123)
##' n <- 16; nout <- 300
##' x <- sort(runif(n))
##' 
##' ##' ## let 'xout' contain n + nout points including nodes 
##' xout <- sort(c(x, runif(nout, min = x[1], max = x[n])))
##' res  <- cardinalBasis_lagrange(x, xout = xout)
##' 
##' matplot(res$x, res$CB, type = "n", main = "Cardinal Basis")
##' abline(v = x, h = 1.0, col = "gray")
##' points(x = x, y = rep(0, n), pch = 21, col = "black", lwd = 2, bg = "white")
##' matlines(res$x, res$CB, type = "l")
##'
##' ## Lebesgue's function is constant = 1.0: check it
##' L <- apply(res$CB, 1, function(x) sum(abs(x)))
##' range(L)

cardinalBasis_lagrange <- function(x, xout) {
    
    res <- cardinalBasis_ceschino(x, xout, cubic = FALSE, deriv = 0) 
    res[["method"]] <- "Lagrange (piecewise linear) interpolation"
    res
}
