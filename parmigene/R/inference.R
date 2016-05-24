# Copyright 2011 Gabriele Sales <gabriele.sales@unipd.it>
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


aracne.a <- function(mi, eps=0.05) {
  if (!is.matrix(mi))
    stop("mi must be a matrix.")

  n <- nrow(mi)
  if (ncol(mi) != n)
    stop("mi must be a square matrix.")

  if (eps <= 0)
    stop("eps must be positive.")

  .aracne(mi, n, eps, 1)
}

aracne.m <- function(mi, tau=0.15) {
  if (!is.matrix(mi))
    stop("mi must be a matrix.")

  n <- nrow(mi)
  if (ncol(mi) != n)
    stop("mi must be a square matrix.")

  if (tau <= 0)
    stop("tau must be positive.")

  .aracne(mi, n, 0, 1-tau)
}

clr <- function(mi) {
  if (!is.matrix(mi))
    stop("mi must be a matrix.")

  n <- nrow(mi)
  if (ncol(mi) != n)
    stop("mi must be a square matrix.")

  .clr(mi, n)
}

mrnet <- function(mi) {
  if (!is.matrix(mi))
    stop("mi must be a matrix.")

  n <- nrow(mi)
  if (ncol(mi) != n)
    stop("mi must be a square matrix.")

  .mrnet(mi, n)
}
