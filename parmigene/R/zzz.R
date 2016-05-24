# Copyright 2010-2012 Gabriele Sales <gabriele.sales@unipd.it>
#
#
# This file is part of parmigene.
#
# knnmi is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License
# version 3 as published by the Free Software Foundation.
#
# knnmi is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public
# License along with parmigene. If not, see <http://www.gnu.org/licenses/>.


.mi_single <- function(x, y, k, noise)
  .C("mi_single",
     as.double(x),
     as.double(y),
     as.integer(length(x)),
     as.integer(k),
     as.double(noise),
     res = double(1),
     PACKAGE = "parmigene",
     DUP = TRUE)$res

.mi_cross <- function(xs, ys, k, noise) {
  h1 <- nrow(xs)
  h2 <- nrow(ys)
  w  <- ncol(xs)
  
  res <- .C("mi_cross",
            as.double(t(xs)),
            as.integer(h1),
            as.double(t(ys)),
            as.integer(h2),
            as.integer(w),
            as.integer(k),
            as.double(noise),
            res = double(h1*h2),
            PACKAGE = "parmigene",
            DUP = TRUE)$res

  m <- t(matrix(res, nrow=h2))
  rownames(m) <- rownames(xs)
  colnames(m) <- rownames(ys)
  m
}

.mi_all <- function(xs, k, noise) {
  h <- nrow(xs)
  w <- ncol(xs)
  
  res <- .C("mi_all",
            as.double(t(xs)),
            as.integer(h),
            as.integer(w),
            as.integer(k),
            as.double(noise),
            res = double(h*h),
            PACKAGE = "parmigene",
            DUP = TRUE)$res

  m <- matrix(res, nrow=h)
  colnames(m) <- rownames(xs)
  rownames(m) <- rownames(xs)
  m
}

.aracne <- function(mis, n, eps, eta) {
  res <- .C("aracne",
            as.double(t(mis)),
            as.integer(n),
            as.double(eps),
            as.double(eta),
            res=double(n*n),
            PACKAGE = "parmigene",
            DUP = FALSE)$res

  m <- matrix(res, nrow=n)
  colnames(m) <- rownames(mis)
  rownames(m) <- rownames(mis)
  m
}

.clr <- function(mis, n) {
  res <- .C("clr",
            as.double(t(mis)),
            as.integer(n),
            res = double(n*n),
            PACKAGE = "parmigene",
            DUP = FALSE)$res

  m <- matrix(res, nrow=n)
  colnames(m) <- rownames(mis)
  rownames(m) <- rownames(mis)
  m
}

.mrnet <- function(mis, n) {
    res <- .C("mrnet",
              as.double(t(mis)),
              as.integer(n),
              res = double(n*n),
              PACKAGE = "parmigene",
              DUP = FALSE)$res

  m <- matrix(res, nrow=n)
  colnames(m) <- rownames(mis)
  rownames(m) <- rownames(mis)
  m
}
