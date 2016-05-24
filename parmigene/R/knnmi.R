# Copyright 2010-2011 Gabriele Sales <gabriele.sales@unipd.it>
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


.check.matrix <- function(x, k, name) {
  if ((!is.matrix(x) && !is.data.frame(x)) || nrow(x) < 2) {
    stop(paste(name, "must be a multi-row matrix or a data.frame"))
  } else if (ncol(x) <= k) {
    stop(paste(name, " has too few columns (must be > ", k, ")", sep=""))
  }
}

knnmi <- function(x, y, k=3, noise=1e-12) {
  if (!is.vector(x))
    stop("x must be a vector.")

  if (!is.vector(y) || length(x) != length(y))
    stop("x is a vector; y should be a vector of the same length")

  if (length(x) <= k) {
    stop(paste("vector length should be greater than", k))
  } else if (k < 2) {
    stop("k must be >= 2.")
  }

  .mi_single(x, y, k, noise)
}

knnmi.cross <- function(mat1, mat2, k=3, noise=1e-12) {
  .check.matrix(mat1, k, "mat1")
  .check.matrix(mat2, k, "mat2")

  if (ncol(mat1) != ncol(mat2)) {
    stop("mat1 and mat2 must have the same number of columns.")
  } else if (k < 2) {
    stop("k must be >= 2.")
  }

  .mi_cross(mat1, mat2, k, noise)
}

knnmi.all <- function(mat, k=3, noise=1e-12) {
  .check.matrix(mat, k, "xs")

  if (k < 2)
    stop("k must be >= 2.")

  .mi_all(mat, k, noise)
}
