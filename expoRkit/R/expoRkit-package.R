##' R-interface to the Fortran package Expokit for matrix
##' exponentiation.
##'
##' \tabular{ll}{
##' Package: \tab expoRkit \cr
##' Type: \tab Package \cr
##' Version: \tab 0.9 \cr
##' Date: \tab 2012-10-08\cr
##' License: \tab GPL (>= 2)\cr
##' LazyLoad: \tab yes \cr
##' }
##'
##' Expokit is an efficient Fortran implementation for computing the
##' matrix exponential, or rather, its action on a vector, for large
##' sparse matrices. This can also be understood as computing the
##' solution of a system of linear ordinary first order differential
##' equations. This package provides an R-interface to some of the
##' Fortran subroutines from Expokit.
##'
##' The Fortran package was developed by Roger B. Sidje, see
##' \url{http://www.maths.uq.edu.au/expokit/}. Niels Richard Hansen
##' adapted the package for use with R and wrote the R
##' interface. Permission to distribute the Expokit source under GPL was
##' obtained from Roger B. Sidje.
##' 
##' @title Expokit in R
##' @name expoRkit
##' @references Sidje, R. B. (1998) Expokit. Software Package for Computing Matrix
##' Exponentials. ACM Trans. Math. Softw. 24(1), 130-156.
##' @seealso \code{\link{expv}}
##' @author Roger B. Sidje \email{roger.b.sidje@@ua.edu}, Niels Richard Hansen
##' \email{Niels.R.Hansen@@math.ku.dk}.
##' @keywords package 
##' @docType package
NULL
