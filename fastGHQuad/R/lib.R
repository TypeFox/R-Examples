#' Evaluate Hermite polynomial at given location
#' 
#' Evaluate Hermite polynomial of given degree at given location. This function
#' is provided for demonstration/teaching purposes; this method is not used by
#' gaussHermiteData. It is numerically unstable for high-degree polynomials.
#' 
#' 
#' @param x Vector of location(s) at which polynomial will be evaluated
#' @param n Degree of Hermite polynomial to compute
#' @return Vector of length(x) values of Hermite polynomial
#' @author Alexander W Blocker \email{ablocker@@gmail.com}
#' @export
#' @seealso \code{\link{gaussHermiteData}}, \code{\link{aghQuad}},
#' \code{\link{ghQuad}}
#' @keywords math
evalHermitePoly <- function(x, n) {
    .Call("evalHermitePoly", x, n, PACKAGE="fastGHQuad")
}



#' Find real parts of roots of polynomial
#' 
#' Finds real parts of polynomial's roots via eigendecomposition of companion
#' matrix. This method is not used by gaussHermiteData. Only the real parts of
#' each root are retained; this can be useful if the polynomial is known a
#' priori to have all roots real.
#' 
#' 
#' @param c Coefficients of polynomial
#' @return Numeric vector containing the real parts of the roots of the
#' polynomial defined by c
#' @author Alexander W Blocker \email{ablocker@@gmail.com}
#' @export
#' @seealso \code{\link{gaussHermiteData}}, \code{\link{aghQuad}},
#' \code{\link{ghQuad}}
#' @keywords math
findPolyRoots <- function(c) {
    .Call("findPolyRoots", c, PACKAGE="fastGHQuad")
}



#' Get coefficient of Hermite polynomial
#' 
#' Calculate coefficients of Hermite polynomial using recursion relation. This
#' function is provided for demonstration/teaching purposes; this method is not
#' used by gaussHermiteData. It is numerically unstable for high-degree
#' polynomials.
#' 
#' 
#' @param n Degree of Hermite polynomial to compute
#' @return Vector of (n+1) coefficients from requested polynomial
#' @author Alexander W Blocker \email{ablocker@@gmail.com}
#' @export
#' @seealso \code{\link{gaussHermiteData}}, \code{\link{aghQuad}},
#' \code{\link{ghQuad}}
#' @keywords math
hermitePolyCoef <- function(n) {
    .Call("hermitePolyCoef", n, PACKAGE="fastGHQuad")
}



#' Compute Gauss-Hermite quadrature rule
#' 
#' Computes Gauss-Hermite quadrature rule of requested order using Golub-Welsch
#' algorithm. Returns result in list consisting of two entries: x, for nodes,
#' and w, for quadrature weights. This is very fast and numerically stable,
#' using the Golub-Welsch algorithm with specialized eigendecomposition
#' (symmetric tridiagonal) LAPACK routines. It can handle quadrature of order
#' 1000+.
#' 
#' This function computes the Gauss-Hermite rule of order n using the
#' Golub-Welsch algorithm. All of the actual computation is performed in C/C++
#' and FORTRAN (via LAPACK). It is numerically-stable and extremely
#' memory-efficient for rules of order 1000+.
#' 
#' @param n Order of Gauss-Hermite rule to compute (number of nodes)
#' @return A list containing: \item{x}{the n node positions for the requested
#' rule} \item{w}{the w quadrature weights for the requested rule}
#' @author Alexander W Blocker \email{ablocker@@gmail.com}
#' @export
#' @seealso \code{\link{aghQuad}}, \code{\link{ghQuad}}
#' @references Golub, G. H. and Welsch, J. H. (1969). Calculation of Gauss
#' Quadrature Rules. Mathematics of Computation 23 (106): 221-230
#' 
#' Liu, Q. and Pierce, D. A. (1994). A Note on Gauss-Hermite Quadrature.
#' Biometrika, 81(3) 624-629.
#' @keywords math
gaussHermiteData <- function(n) {
    .Call("gaussHermiteData", n, PACKAGE="fastGHQuad")
}
