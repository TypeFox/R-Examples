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


##
Pconstruct <- function(theta){
  n <- length(theta)
  d <- (1 + sqrt(1 + 8 * n)) / 2
  A <- matrix(0, nrow = d, ncol = d)
  A[lower.tri(A)] <- theta
  diag(A) <- 1
  Q <- A %*% t(A)
  P <- cov2cor(Q)
  P
}
##
Pdeconstruct <- function(P){
  A <- t(chol(P))
  Adiag <- diag(diag(A))
  Astar <- solve(Adiag) %*% A
  Astar[lower.tri(Astar)]
}
## Empirical (cumulative) distribution function
edf <- function(v, adjust = FALSE){
  original <- v
  v <- sort(v)
  vv <- cumsum(!duplicated(v))
  repeats <- tapply(v, v, length)
  add <- rep(cumsum(repeats - 1.), repeats)
  df <- (vv + add)/(length(vv) + as.numeric(adjust))
  as.numeric(df[rank(original)])
}
##
ESnorm <- function(p, mu = 0, sd = 1){
  mu + sd * dnorm(qnorm(p)) / (1 - p)
}
##
ESst <- function(p, mu = 0, sd = 1, df, scale = FALSE){
  ES <- (dt(qt(p, df), df)/(1 - p)) * ((df + (qt(p, df))^2) / (df - 1))
  if (scale) ES <- ES * sqrt((df - 2) / df)
  mu + sd * ES
}
