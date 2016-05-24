#   R package for Singular Spectrum Analysis
#   Copyright (c) 2009 Anton Korobeynikov <asl@math.spbu.ru>
#   'svd' R package
#   Copyright (c) 2015 Anton Korobeynikov <asl@math.spbu.ru>
#
#   This program is free software; you can redistribute it
#   and/or modify it under the terms of the GNU General Public
#   License as published by the Free Software Foundation;
#   either version 2 of the License, or (at your option)
#   any later version.
#
#   This program is distributed in the hope that it will be
#   useful, but WITHOUT ANY WARRANTY; without even the implied
#   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#   PURPOSE.  See the GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public
#   License along with this program; if not, write to the
#   Free Software Foundation, Inc., 675 Mass Ave, Cambridge,
#   MA 02139, USA.

#   Routines for external matrix stuff

extmat.ncol <- function(X)
  .extmat.ncol(X@.xData)

extmat.nrow <- function(X)
  .extmat.nrow(X@.xData)

.extmat.ncol <- function(X) {
  .Call("extmat_ncol", X)
}

.extmat.nrow <- function(X) {
  .Call("extmat_nrow", X)
}

is.extmat <- function(X)
  is(X, "extmat") && .is.extmat(X@.xData)

.is.extmat <- function(X) {
  .Call("is_extmat", X)
}

ematmul <- function(emat, v, transposed = FALSE) {
  storage.mode(v) <- "double";
  storage.mode(transposed) <- "logical";
  .Call("ematmul", emat@.xData, v, transposed);
}

.ematmul <- function(emat, v, transposed = FALSE) {
  storage.mode(v) <- "double";
  storage.mode(transposed) <- "logical";
  .Call("ematmul_unchecked", emat, v, transposed);
}

extmat <- function(mul, tmul, nrow, ncol,
                   env = parent.frame()) {
  new("extmat",
      .Call("initialize_rextmat",
            match.fun(mul), match.fun(tmul),
            as.integer(nrow), as.integer(ncol),
            env))
}

# S4 weirdness
setClass("extmat", contains = "externalptr")
  
setMethod("as.matrix", "extmat",
          function(x) as(x, "matrix"))
setMethod("as.array",  "extmat", function(x) as(x, "matrix"))
as.array.extmat <- as.matrix.extmat <- function(x, ...) as(x, "matrix")

setMethod("as.vector", "extmat",
          function(x, mode) as.vector(as(x, "matrix"), mode))
as.vector.extmat <- function(x, mode) as.vector(as(x, "matrix"), mode)

setMethod("as.numeric", "extmat",
          function(x, ...) as.numeric(as.vector(x)))
setMethod("as.integer", "extmat",
          function(x, ...) as.integer(as.vector(x)))
setMethod("as.logical", "extmat",
          function(x, ...) as.logical(as.vector(x)))

setAs("extmat", "matrix", function(from) from %*% diag(nrow = ncol(from)))

setMethod("dim", "extmat",
          function(x) c(.extmat.nrow(x@.xData), .extmat.ncol(x@.xData)), valueClass = "integer")
setMethod("length", "extmat", function(x) prod(dim(x)))

setMethod("t", "extmat",
          function(x) {
            emat <- x@.xData

            extmat(mul = function(v) .ematmul(emat, v, transposed = TRUE),
                   tmul = function(v) .ematmul(emat, v, transposed = FALSE),
                   nrow = ncol(x), ncol = nrow(x))
          })

setMethod("%*%", signature(x = "extmat", y = "numeric"),
          function(x, y) {
            dim(y) <-
              if (ncol(x) == (n <- length(y))) c(n, 1L) else c(1L, n)
            x %*% y
          })
setMethod("%*%", signature(x = "numeric", y = "extmat"),
          function(x, y) {
            dim(x) <-
              if (nrow(y) == (n <- length(x))) c(1L, n) else c(n, 1L)
            x %*% y
          })

setMethod("%*%", signature(x = "extmat", y = "matrix"),
          function(x, y) {
            if (nrow(y) != ncol(x))
              stop("non-conformable arguments")
            res <- apply(y, 2, .ematmul, emat = x@.xData, transposed = FALSE)
            dim(res) <- c(nrow(x), ncol(y))
            res
          })
setMethod("%*%", signature(x = "matrix", y = "extmat"),
          function(x, y) {
            if (nrow(y) != ncol(x))
              stop("non-conformable arguments")
            res <- t(apply(x, 1, .ematmul, emat = y@.xData, transposed = TRUE))
            dim(res) <- c(nrow(x), ncol(y))
            res
          })

setMethod("%*%", signature(x = "extmat", y = "extmat"),
          function(x, y) {
            res <- apply(diag(ncol(y)), 2,
                         function(u, emat.x, emat.y)
                           .ematmul(emat = emat.x, .ematmul(emat = emat.y, u)),
                         emat.x = x@.xData, emat.y = y@.xData)
            dim(res) <- c(nrow(x), ncol(y))
            res
          })

setMethod("crossprod", signature(x = "extmat", y = "extmat"),
          function(x, y) {
            t(x) %*% y
          })
setMethod("crossprod", signature(x = "extmat", y = "ANY"),
          function(x, y) {
            t(x) %*% y
          })
setMethod("crossprod", signature(x = "ANY", y = "extmat"),
          function(x, y) {
            t(t(y) %*% x)
          })
setMethod("crossprod", signature(x = "extmat", y = "missing"),
          function(x, y) {
            crossprod(x, x)
          })

setMethod("tcrossprod", signature(x = "extmat", y = "extmat"),
          function(x, y) {
            t(y %*% t(x))
          })
setMethod("tcrossprod", signature(x = "extmat", y = "ANY"),
          function(x, y) {
            t(y %*% t(x))
          })
setMethod("tcrossprod", signature(x = "ANY", y = "extmat"),
          function(x, y) {
            x %*% t(y)
          })
setMethod("tcrossprod", signature(x = "extmat", y = "missing"),
          function(x, y) {
            tcrossprod(x, x)
          })
