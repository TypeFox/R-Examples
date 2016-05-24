#################################################################################
##
##   R package Copula by Jun Yan and Ivan Kojadinovic Copyright (C) 2009
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

###################################################################
## EV test based on Cn (local power)
###################################################################

localpowerC <- function(copC, copD, N=1000, n=1000, step=5, ndelta=4, alpha = 0.05)
{
  p <- 2
  x <- rCopula(n, copC)
  u <- apply(x,2,rank)/(n+1) ## make pseudo-observations

  ## set r according to recommendations
  r <- 3:5
  nr <- length(r)

  ## grid = pseudo-observations
  g <- u
  m <- n

  ## nonparam
  ##offset <- offsetect <- -1
  ##s0np <- .C("evtest",
  ##           as.double(u),
  ##           as.integer(n),
  ##           as.integer(p),
  ##           as.double(g),
  ##           as.integer(m),
  ##           as.integer(N),
  ##           as.double(1/r),
  ##           as.integer(nr),
  ##           as.double(offset),
  ##           as.double(offsetect),
  ##           s0 = double(N * nr),
  ##           PACKAGE="copula")$s0

  ##s0np <- matrix(s0np, ncol = nr, byrow = TRUE)
  ##s0np <- apply(s0np,1,sum)

  ## param
  der <- dCdu(copC,g)
  dert <- numeric(0)
  for (i in r)
    dert <- c(dert,as.double(dCdu(copC,g^(1/i))))

  pcopC <- pCopula(copC,g)
  pcopD <- pCopula(g, copD)
  pcopCt <- numeric(0)
  delta.term <- numeric(0)
  for (i in r)
    {
      pcopCr <- pCopula(g^(1/i), copC)
      pcopCt <- c(pcopCt,pcopCr)
      delta.term <- c(delta.term,
                      i *  pcopCr^(i - 1) * (pCopula(g^(1/i), copD) -  pcopCr)
                      - (pcopD - pcopC))
    }

  s0 <- .C("evtest_LP",
           as.double(x),
           as.integer(n),
           as.integer(p),
           as.double(g),
           as.integer(m),
           as.integer(N),
           as.double(1/r),
           as.integer(nr),
           s0 = double(N * ndelta),
           as.double(pcopCt),
           as.double(der),
           as.double(dert),
           as.double(delta.term),
           as.double(step),
           as.double(ndelta),
           PACKAGE="copula")$s0

  s0 <- matrix(s0, ncol = ndelta)

  q <- sort(s0[,1])[(1- alpha) * N] #critical value

  localpow <- numeric(ndelta - 1)
  for (i in 2:ndelta)
    localpow[i-1] <- mean(s0[,i] >= q)

  #return(list(s0np=s0np,s0=s0,lp=c(alpha,localpow)))
  return(c(alpha,localpow))
}

###################################################################
## EV test based on An (local power)
###################################################################

AD <- function(copD,y)
  {
    integrand <- function(x) {(pCopula(cbind(x^(1-y),x^y), copD) - (x > exp(-1)))/ (x * log(x))}
    exp(-0.57721566 + integrate(integrand, 0,1)$value)
  }

localpowerA <- function(copC, copD, N=1000, n=1000, step=5, ndelta=4, alpha = 0.05)
{
  x <- rCopula(n, copC)
  u <- apply(x,2,rank)/(n+1) ## make pseudo-observations

  #s0np <- .C("evtestA",
  #           as.double(u[,1]),
  #           as.double(u[,2]),
  #           as.integer(n),
  #           as.double(u[,1]),
  #           as.double(u[,2]),
  #           as.integer(n),
  #           as.integer(1),
  #           as.integer(N),
  #           s0 = double(N),
  #           PACKAGE="copula")$s0

  der <- dCdu(copC,u)
  loguv <- log(u[,1]*u[,2])
  g <- log(u[,2])/loguv #grid (size n)
  AC <- A(copC,g)
  ADg <- numeric(n)
  for (i in 1:n)
    ADg[i] <- AD(copD,g[i])
  delta.term <- pCopula(u, copD) - pCopula(u, copC) -
    exp(loguv * AC) * loguv * AC * (log(ADg) - log(AC))

  s0 <- .C("evtestA_LP",
           as.double(u[,1]),
           as.double(u[,2]),
           as.integer(n),
           as.double(u[,1]),
           as.double(u[,2]),
           as.integer(n),
           as.integer(N),
           s0 = double(N * ndelta),
           as.double(AC),
           as.double(der[,1]),
           as.double(der[,2]),
           as.double(delta.term),
           as.double(step),
           as.double(ndelta),
           PACKAGE="copula")$s0

  s0 <- matrix(s0, ncol = ndelta)

  q <- sort(s0[,1])[(1- alpha) * N] #critical value

  localpow <- numeric(ndelta - 1)
  for (i in 2:ndelta)
    localpow[i-1] <- mean(s0[,i] >= q)

  #return(list(s0np=s0np,s0=s0,lp=c(alpha,localpow)))
  return(c(alpha,localpow))
}
