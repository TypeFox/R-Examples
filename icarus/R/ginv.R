# This function is directly taken from MASS package so as to avoid
# having too many dependancies (code on GPLv3 license, I do not hold
# copyright on this function).
# copyright (C) 1994-2014 W. N. Venables and B. D. Ripley
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

ginv <- function(X, tol = sqrt(.Machine$double.eps))
{
  #
  # based on suggestions of R. M. Heiberger, T. M. Hesterberg and WNV
  #
  if(length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X)))
    stop("'X' must be a numeric or complex matrix")
  if(!is.matrix(X)) X <- as.matrix(X)
  Xsvd <- svd(X)
  if(is.complex(X)) Xsvd$u <- Conj(Xsvd$u)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
  if (all(Positive)) Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
  else if(!any(Positive)) array(0, dim(X)[2L:1L])
  else Xsvd$v[, Positive, drop=FALSE] %*% ((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop=FALSE]))
}
