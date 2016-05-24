# Copyright (C) 2014
# This function was borrowed from the package RGCCA tau.estimate.R with small changes

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


tau.estimate <- function(x) {
  if (is.matrix(x) == TRUE && is.numeric(x) == FALSE) 
    stop("Data matrices in the list must be numeric")
  p <- NCOL(x)
  n <- NROW(x)
  covm <- cov(x)
  corm <- cor(x)
  
  # the matrices are centered and scaled on columns to compute the crossprods below
  xs <- scale(x, center = TRUE, scale = TRUE)
  
  # we use the Strimmer formula
  v <- (n/((n - 1)^3)) * (crossprod(xs^2) - 1/n * (crossprod(xs))^2)
  diag(v) <- 0
  m <- matrix(rep(apply(xs^2, 2, mean), p), p, p)
  I <- diag(NCOL(x))
  d <- (corm - I)^2
  
  tau <- (sum(v))/sum(d)
  tau <- max(min(tau, 1), 0)
  return(tau)
}
