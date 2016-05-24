################################################################################
################################################################################
## Generate data that fit a given direct preferences matrix
## by particle swarm optimization

# Copyright 2014 Julien Boelaert.
# 
# This file is part of revealedPrefs.
# 
# revealedPrefs is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# revealedPrefs is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with revealedPrefs.  If not, see <http://www.gnu.org/licenses/>.

## Cost function used in PSO
optCost <- function(x, objective.mat, afriat.par) {
  matx <- matrix(x[1:(length(x)/2)], nrow= nrow(objective.mat))
  matp <- matrix(x[((length(x)/2)+1):length(x)], nrow= nrow(objective.mat))  
  px <- matp %*% t(matx)
  px <- afriat.par * diag(px) - px
  
  res <- matrix(NA, nrow= nrow(objective.mat), ncol= ncol(objective.mat))
  res[objective.mat == 1] <- abs(px[objective.mat == 1])
  res[objective.mat == 2] <- max(0, 1e-6 - px[objective.mat == 2])
  res[objective.mat == 0] <- max(0, 1e-6 + px[objective.mat == 0])
  sum(res)
}

## function to generate simulated data that fit direct preference matrix obj.mat
## with n.goods types of goods
simPrefs <- function(pref.mat, ngoods, afriat.par= 1, 
                     qmin= 0, qmax= 1, pmin= 0, pmax= 1,
                     maxit= 100, verbose= FALSE) {
  library(pso)
  if (length(afriat.par) > 1 | afriat.par > 1 | afriat.par < 0)
    stop("'afriat.par' must be a real value between 0 and 1.\n")
  res <- pso::psoptim(par = rep(NA, 2 * nrow(pref.mat) * ngoods), 
                      fn = function(x) optCost(x, pref.mat, afriat.par),
                      lower = c(rep(qmin, ngoods), rep(pmin, ngoods)), 
                      upper = c(rep(qmax, ngoods), rep(pmax, ngoods)),
                      control= list(trace= ifelse(verbose, 1, 0), maxit= maxit))
  if (res$value == 0)
    return(list(x= matrix(res$par[1:(length(res$par)/2)], nrow= nrow(pref.mat)), 
                p= matrix(res$par[((length(res$par)/2)+1):length(res$par)], 
                          nrow= nrow(pref.mat))))
  warning("No solution found. Try again, or with more iterations, or with another matrix.")
  return(list(x= NULL, p= NULL))
}

################################################################################
## Iteratively generate random GARP-consistent observations
simGarp <- function(nobs, ngoods, afriat.par= 1, maxit= 10 * nobs, 
                    qmin= 0, qmax= 1, pmin= 0, pmax= 1) {
  if (length(afriat.par) > 1 | afriat.par > 1 | afriat.par < 0)
    stop("'afriat.par' must be a real value between 0 and 1.\n")
  
  res <- .Call("SimAxiom", nobs, ngoods, afriat.par, maxit, 
               pmin, pmax, qmin, qmax, "GARP", PACKAGE= "revealedPrefs")
  if (res$nobs < nobs & res$iter == maxit)
    warning("Max iterations reached before completing requested size.")
  res
}

################################################################################
## Iteratively generate random SARP-consistent observations
simSarp <- function(nobs, ngoods, afriat.par= 1, maxit= 10 * nobs, 
                    qmin= 0, qmax= 1, pmin= 0, pmax= 1) {
  if (length(afriat.par) > 1 | afriat.par > 1 | afriat.par < 0)
    stop("'afriat.par' must be a real value between 0 and 1.\n")
  
  res <- .Call("SimAxiom", nobs, ngoods, afriat.par, maxit, 
               pmin, pmax, qmin, qmax, "SARP", PACKAGE= "revealedPrefs")
  if (res$nobs < nobs & res$iter == maxit)
    warning("Max iterations reached before completing requested size.")
  res
}

################################################################################
## Iteratively generate random GARP-consistent observations
simWarp <- function(nobs, ngoods, afriat.par= 1, maxit= 10 * nobs, 
                    qmin= 0, qmax= 1, pmin= 0, pmax= 1) {
  if (length(afriat.par) > 1 | afriat.par > 1 | afriat.par < 0)
    stop("'afriat.par' must be a real value between 0 and 1.\n")
  
  res <- .Call("SimAxiom", nobs, ngoods, afriat.par, maxit, 
               pmin, pmax, qmin, qmax, "WARP", PACKAGE= "revealedPrefs")
  if (res$nobs < nobs & res$iter == maxit)
    warning("Max iterations reached before completing requested size.")
  res
}
