## Copyright (C) 2013 Marius Hofert, Bernhard Pfaff
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


## Function for computing equal correlations
equicorr <- function(d, rho){
  d <- as.integer(d)
  if(abs(rho) > 1)
    stop("\nCorrelation is greater than unity in absolute value.\n")
  if(rho < ( - (d - 1.0)^(-1.0)))
    stop(paste("\n'rho' must be at least",  - (d - 1.0)^(-1.0), ".\n"))
  out <- matrix(rho, nrow = d, ncol = d)
  diag(out) <- 1
  out
}

## TODO: deprecate -- use Matrix' nearPD

## Make a matrix positive definite
eigenmeth <- function(mat, delta = 0.001){
  decomp <- eigen(mat)
  Lambda <- decomp$values
  Lambda[Lambda < 0] <- delta
  Gamma <- decomp$vectors
  newmat <- Gamma %*% diag(Lambda) %*% t(Gamma)
  D <- 1/sqrt(diag(newmat))
  diag(D) %*% newmat %*% diag(D)
}

## TODO: should be deprecated -- not useful!

## Spearman rank correlations
Spearman <- function(data, ...){
  out <- cor(data, method = "spearman", ...)
  out
}
## Kendall rank correlations
Kendall <- function(data, ...){
  out <- cor(data, method = "kendall", ...)
  out
}
