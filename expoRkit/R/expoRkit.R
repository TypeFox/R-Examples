### R interface for calling Expokit Fortran routines for computing the
### matrix exponential. This file contains the low level interface to
### the Fortran routines
###
###     Copyright (C) 2012 Niels Richard Hansen.
###
### This program is free software; you can redistribute it and/or modify it
### under the terms of the GNU General Public License as published by the
### Free Software Foundation; either version 2, or (at your option) any
### later version.
###
### These functions are distributed in the hope that they will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
### GNU General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program; if not, a copy is available at
### http://www.r-project.org/Licenses/


##' R wrapper of the Expokit Fortran subroutines __EXPV and __PHIV for
##' sparse matrix exponentiation. In general, these routines compute
##' the solution at time point \eqn{t} of the ODE \deqn{w'(t) = x w(t)
##' + u} with initial condition \eqn{w(0) = v}. 
##' 
##'
##' The \code{Rexpv} function is the low level wrapper of the Fortran
##' subroutines in the Expokit package. It is not intended to be used
##' directly but rather via the \code{\link{expv}} methods. In the
##' call the correct storage format in terms of the vectors \code{a},
##' \code{ia} and \code{ja} has to be specified via the \code{storage}
##' argument. For CCS, \code{ia} contains 1-based row numbers of
##' non-zero elements and \code{ja} contains 1-based pointers to the
##' initial element for each column. For CRS, \code{ja} contains
##' 1-based column numbers of non-zero elements and \code{ia} are
##' 1-based pointers to the initial element for each row. For COO,
##' \code{ia} and \code{ja} contain the 1-based column and row
##' numbers, respectively, for the non-zero elements.
##' @title Expokit __EXPV and __PHIV wrapper.
##' @param a \code{numeric} or \code{complex} non-zero entries in the \eqn{x}-matrix.
##' @param ia \code{integer} index/pointer. Precise meaning depends on
##' storage format.
##' @param ja \code{integer} index/pointer. Precise meaning depends on
##' storage format.
##' @param n dimension of the (square) matrix.
##' @param v \code{numeric} or \code{complex} vector.
##' @param t time. Default \code{1}.
##' @param storage \code{character}, one of \code{'CCS'} (Compressed
##' Column Storage), \code{'CRS'} (Compressed Row Storage) or
##' \code{'COO'} (COOrdinate list). Default \code{'CCS'}.
##' @param u \code{numeric} or \code{complex} vector. Default \code{NULL}.
##' @param anorm A norm of the matrix. Default is the sup-norm.
##' @param Markov \code{logical}, if \code{TRUE} the (transposed)
##' matrix is taken to be an intensity matrix and steps are taken to
##' ensure that the computed result is a probability vector. Default \code{FALSE}. 
##' @param m \code{integer}, the maximum size for the Krylov basis. 
##' @param tol \code{numeric}. A value of 0 (default) means square
##' root of machine eps.
##' @param itrace \code{integer}, 0 (default) means no trace
##' information from Expokit, 1 means print 'happy breakdown', and 2
##' means all trace information printed from Expokit.
##' @param mxstep \code{integer}. Maximum allowable number of
##' integration steps. The value 0 means an infinite number of steps. Default 10000.
##' @return The solution, \eqn{w}, of the ODE as a \code{numeric} or
##' \code{complex} vector of length \eqn{n}.
##' @references Sidje, R. B. (1998) Expokit. Software Package for Computing Matrix
##' Exponentials. ACM Trans. Math. Softw. 24(1), 130-156.
##' @seealso \code{\link[Matrix]{expm}}, \code{\link[expm]{expm}},
##' \code{\link{expv}}
##' @author Niels Richard Hansen \email{Niels.R.Hansen@@math.ku.dk}
##' @examples
##' ### A CCS 4 by 4 real matrix. The last element in 'ja' is the number of
##' ### non-zero elements + 1. 
##' a <- c(-1, 1, -2, -3, 1, 2, -1)
##' ia <- c(1, 3, 2, 4, 1, 2, 3) 
##' ja <- c(1, 3, 5, 6, 8)
##' 
##' v <- c(1, 1, 1, 1)
##' wCCS <- expoRkit:::Rexpv(a, ia, ja, 4, v = v)
##' 
##' ### COO storage instead.
##' ja <- c(1, 1, 2, 2, 3, 4, 4)  
##' wCOO <- expoRkit:::Rexpv(a, ia, ja, 4, v = v, storage = 'COO')
##' 
##' ### CRS storage instead.
##' a <- c(-1, 1, -2, 2, 1, -1, -3)
##' ja <- c(1, 3, 2, 4, 1, 4, 2)
##' ia <- c(1, 3, 5, 7, 8)
##' wCRS <- expoRkit:::Rexpv(a, ia, ja, 4, v = v, storage = 'CRS')
##' 
##' cbind(wCCS, wCOO, wCRS)
##' 
##' stopifnot(all.equal(wCCS, wCOO),
##'           all.equal(wCCS, wCRS),
##'           all.equal(wCRS, wCOO))
##'
##' ## Complex version
##' a <- c(-1, 1i, -2, -3i, 1, 2i, -1)
##' ia <- c(1, 3, 2, 4, 1, 2, 3) 
##' ja <- c(1, 3, 5, 6, 8)
##' 
##' v <- c(1, 1, 1, 1)
##' wCCS <- expoRkit:::Rexpv(a, ia, ja, 4, v = v)
##' 
##' ### COO storage instead.
##' ja <- c(1, 1, 2, 2, 3, 4, 4)  
##' wCOO <- expoRkit:::Rexpv(a, ia, ja, 4, v = v, storage = 'COO')
##' 
##' ### CRS storage instead.
##' a <- c(-1, 1, -2, 2i, 1i, -1, -3i)
##' ja <- c(1, 3, 2, 4, 1, 4, 2)
##' ia <- c(1, 3, 5, 7, 8)
##' wCRS <- expoRkit:::Rexpv(a, ia, ja, 4, v = v, storage = 'CRS')
##' 
##' cbind(wCCS, wCOO, wCRS)
##' 
##' stopifnot(all.equal(wCCS, wCOO),
##'           all.equal(wCCS, wCRS),
##'           all.equal(wCRS, wCOO))
##'
##' @useDynLib expoRkit
Rexpv <- function(a, ia, ja, n, v, t = 1.0, storage = 'CCS', u = NULL,
                  anorm = max(abs(a)), Markov = FALSE, m = 30L, tol =
                  0.0, itrace = 0L, mxstep = 10000L) {
  if (n <= 1)
    stop("The matrix dimension, argument 'n', cannot be 1.")
  m <- as.integer(min(m, n-1))
  sflag <- switch(storage,
                  "CCS" = 1,
                  "CRS" = 2,
                  "COO" = 3,
                  NULL)
  if (is.null(sflag))
    stop("Storage format wrong.")
  if (is.null(u)) {
    if (length(v) != n)
      stop("The length of 'v' does not match the dimension of the matrix.")
    if (is.numeric(a) & !Markov) {
      outlist <- .Fortran("R_DGEXPV",
                          a = as.numeric(a),
                          ia = as.integer(ia),
                          ja = as.integer(ja),
                          n = as.integer(n)[1],
                          nz = as.integer(length(a)),
                          m = m, 
                          t = as.numeric(t)[1],
                          v = as.numeric(v),
                          w = numeric(n),
                          tol = as.numeric(tol)[1],
                          mxstep = as.integer(mxstep),
                          anorm = as.numeric(anorm)[1],
                          itrace = as.integer(itrace)[1],
                          iflag = integer(1),
                          sflag = as.integer(sflag)[1]
                          )
    } else if (is.numeric(a) & Markov) {
      outlist <- .Fortran("R_DMEXPV",
                          a = as.numeric(a),
                          ia = as.integer(ia),
                          ja = as.integer(ja),
                          n = as.integer(n)[1],
                          nz = as.integer(length(a)),
                          m = m, 
                          t = as.numeric(t)[1],
                          v = as.numeric(v),
                          w = numeric(n),
                          tol = as.numeric(tol)[1],
                          mxstep = as.integer(mxstep),
                          anorm = as.numeric(anorm)[1],
                          itrace = as.integer(itrace)[1],
                          iflag = integer(1),
                          sflag = as.integer(sflag)[1]
                          )
    } else if (is.complex(a) & !Markov) {
      outlist <- .Fortran("R_ZGEXPV",
                          a = as.complex(a),
                          ia = as.integer(ia),
                          ja = as.integer(ja),
                          n = as.integer(n)[1],
                          nz = as.integer(length(a)),
                          m = m, 
                          t = as.numeric(t)[1],
                          v = as.complex(v),
                          w = complex(n),
                          tol = as.numeric(tol)[1],
                          mxstep = as.integer(mxstep),
                          anorm = as.numeric(anorm)[1],
                          itrace = as.integer(itrace)[1],
                          iflag = integer(1),
                          sflag = as.integer(sflag)[1]
                          )
    } else {
      stop("Markov condition cannot be enforced for a complex matrix.")
    }
  } else {
    if (length(v) != n)
      stop("The length of 'v' does not match the dimension of the matrix.")
    if (length(u) != n)
      stop("The length of 'u' does not match the dimension of the matrix.")
    if (is.numeric(a)) {
      outlist <- .Fortran("R_DGPHIV",
                          a = as.numeric(a),
                          ia = as.integer(ia),
                          ja = as.integer(ja),
                          n = as.integer(n)[1],
                          nz = as.integer(length(a)),
                          m = m, 
                          t = as.numeric(t)[1],
                          u = as.numeric(u),
                          v = as.numeric(v),
                          w = numeric(n),
                          tol = as.numeric(tol)[1],
                          mxstep = as.integer(mxstep),
                          anorm = as.numeric(anorm)[1],
                          itrace = as.integer(itrace)[1],
                          iflag = integer(1),
                          sflag = as.integer(sflag)[1]
                          )
    } else if (is.complex(a)) {
      outlist <- .Fortran("R_ZGPHIV",
                          a = as.complex(a),
                          ia = as.integer(ia),
                          ja = as.integer(ja),
                          n = as.integer(n)[1],
                          nz = as.integer(length(a)),
                          m = m, 
                          t = as.numeric(t)[1],
                          u = as.complex(u),
                          v = as.complex(v),
                          w = complex(n),
                          tol = as.numeric(tol)[1],
                          mxstep = as.integer(mxstep),
                          anorm = as.numeric(anorm)[1],
                          itrace = as.integer(itrace)[1],
                          iflag = integer(1),
                          sflag = as.integer(sflag)[1]
                          )
    } else {
      stop("The matrix entries must be 'numeric' or 'complex'.")
    }
  }
    if (outlist$iflag != 0) 
      stop(paste("Problem. Exit flag", outlist$iflag, "from Expokit"))
    outlist$w
}

##' R wrapper of the Expokit subroutines __PADM for dense matrix
##' exponentiation via the Pade approximation.
##'
##' The underlying Fortran routines compute the matrix exponential
##' using a combination of scaling and squaring and the irreducible
##' rational Pade approximation.
##' @title Pade approximation of dense matrix exponential. 
##' @param x \code{numeric} or \code{complex} matrix. 
##' @param t time. Default 1. 
##' @param order \code{integer}, the order of the Pade
##' approximation. The default (6) is usually enough.
##' @return The matrix exponential, \eqn{\exp(tx)}{exp(tx)}, as a
##' \code{numeric} or \code{complex} matrix.
##' @references Sidje, R. B. (1998) Expokit. Software Package for Computing Matrix
##' Exponentials. ACM Trans. Math. Softw. 24(1), 130-156.
##' @seealso \code{\link[Matrix]{expm}}, \code{\link[expm]{expm}},
##' \code{\link{expv}}
##' @author Niels Richard Hansen \email{Niels.R.Hansen@@math.ku.dk}
##' @export
padm <- function(x, t = 1.0, order = 6L) {
  x <- as.matrix(x)
  n <- as.integer(nrow(x))
  if (n != ncol(x))
    stop("Matrix argument 'x' must be square.")
  order <- as.integer(order)[1]
  lwsp <- 4L*n*n + order + 1L
  if (is.numeric(x)) {
    outlist <- .Fortran("DGPADM",
                        ideg = order,
                        m = n,
                        t = as.numeric(t)[1],
                        H = as.numeric(x),
                        ldh = n,  ## Any cases where this is not n?
                        wsp = numeric(lwsp),
                        lwsp = lwsp,
                        ipiv = integer(n),
                        iexph = integer(1),
                        ns = integer(1),
                        iflag = integer(1)
                        )
  } else if (is.complex(x)) {
    outlist <- .Fortran("ZGPADM",
                        ideg = order,
                        m = n,
                        t = as.numeric(t)[1],
                        H = as.complex(x),
                        ldh = n,  ## Any cases where this is not n?
                        wsp = complex(lwsp),
                        lwsp = lwsp,
                        ipiv = integer(n),
                        iexph = integer(1),
                        ns = integer(1),
                        iflag = integer(1)
                        )
  } else {
    stop("Argument 'x' must be a 'numeric' or 'complex' matrix.")
  }
  if (outlist$iflag != 0) 
    stop(paste("Problem. Exit flag", outlist$iflag, "from Expokit."))
  matrix(outlist$wsp[seq(outlist$iexph, outlist$iexph + n*n - 1)], n, n)   
}
