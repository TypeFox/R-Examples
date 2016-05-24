#################################################################################
##
##   R package Copula by Jun Yan and Ivan Kojadinovic Copyright (C) 2008
##
##   This file is part of the R package copula.
##
##   The R package copula is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package copula is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package copula. If not, see <http://www.gnu.org/licenses/>.
##
#################################################################################

#################################################################################

## variance/covariance of the pseudo-likelihood estimator

influ.terms <- function(u, influ, q)
{
  p <- ncol(u)
  n <- nrow(u)

  o <- ob <- matrix(0,n,p)
  for (i in 1:p)
    {
      o[,i] <- order(u[,i], decreasing=TRUE)
      ob[,i] <- rank(u[,i])
    }

  out <- matrix(0,n,q)
  for (i in 1:p)
      out <- out + rbind(rep(0,q),apply(influ[[i]][o[,i],,drop=FALSE],2,cumsum))[n + 1 - ob[,i],,drop=FALSE] / n
  return(out)
}

## cop is the FITTED copula
## u are the available pseudo-observations

varPL <- function(cop,u)
  {
    p <- cop@dimension
    n <- nrow(u)

    ## influence: second part
    ## integrals computed from the original pseudo-obs u by Monte Carlo

    dcop <- dcopwrap(cop,u) ## wrapper
    influ0 <- scale(derPdfWrtParams(cop,u)/dcop,scale=FALSE)
    derArg <- derPdfWrtArgs(cop,u)/dcop

    influ <- vector("list",p)
    for (i in 1:p)
        influ[[i]] <- influ0 * derArg[,i]

    ## expectation
    q <- length(cop@parameters)
    e <- crossprod(influ0)
    e <- e/n

    return(var((influ0 - influ.terms(u,influ,q)) %*% solve(e)))
  }


#################################################################################

## variance of the estimator based on Kendall's tau

## cop is the FITTED copula
## u are the available pseudo-observations

varKendall <- function(cop,u)
{
  return(var(4 * (2 * pCopula(u, cop) - u[,1] - u[,2]) / dTau(cop)))
}

varKendallNP <- function(cop,u)
{
  n <- nrow(u)
  ec <- numeric(n)
  for (i in 1:n)
    ec[i] <- sum(u[,1] <= u[i,1] & u[,2] <= u[i,2])/n

  return(var(4 * (2 * ec - u[,1] - u[,2]) / dTau(cop)))
}

##############################################################################

## variance of the estimator based on Spearman's rho

## additional influence terms

add.terms <- function(u)
{
  n <- nrow(u)
  o1 <- order(u[,1], decreasing=TRUE)
  o1b <- rank(u[,1])
  o2 <- order(u[,2], decreasing=TRUE)
  o2b <- rank(u[,2])
  return(c(0,cumsum(u[o1,2]))[n + 1 - o1b] / n  +
         c(0,cumsum(u[o2,1]))[n + 1 - o2b] / n )
}

## cop is the FITTED copula
## u are the available pseudo-observations

varSpearman <- function(cop,u)
  {
    return(var((12 * (u[,1] * u[,2] + add.terms(u))) / dRho(cop)))
  }

##############################################################################
