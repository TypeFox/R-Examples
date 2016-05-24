#******************************************************************************* 
#
# Estimation for Multivariate Normal Data with Monotone Missingness
# Copyright (C) 2007, University of Cambridge
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
# Questions? Contact Robert B. Gramacy (bobby@statslab.cam.ac.uk)
#
#*******************************************************************************


## kl.norm
##
## Kulback-Leibler divergence between two multivariate normal
## distributions specified by their mean vectors and
## covariance matrices -- if a covariance matrix is not
## positive definite, then the nearst posdef one is found and
## used in the distance calculation

`kl.norm` <-
function(mu1, S1, mu2, S2, quiet=FALSE, symm=FALSE)
{
  N <- length(mu1)

  ## check that the mean vectors have the same length
  if(length(mu2) != N) stop("must have length(mu1) == length(mu2)")

  ## check the covar matrices have same dims as the mean
  if(ncol(S1) != N || nrow(S1) != N)
    stop("must have nrow(S1) == ncol(S1) == length(mu1)")
  if(ncol(S2) != N || nrow(S2) != N)
    stop("must have nrow(S2) == ncol(S2) == length(mu2)")

  ## cholesky
  S1c <- chol(S1)
  S2c <- chol(S2)

  ##
  ## distance calculation in parts
  ##

  ## calculate the determinants in log space
  ## ld2 <- determinant(S2, logarithm=TRUE)
  ## if(ld2$sign == -1) { warning("S2 is not posdef"); return(Inf) }
  ld2 <- 2*sum(log(diag(S2c)))
  ##ld1 <- determinant(S1, logarithm=TRUE)
  ## if(ld1$sign == -1) { warning("S1 is not posdef"); return(Inf) }
  ld1 <- 2*sum(log(diag(S1c)))
  ##ldet <- ld2$modulus[1] - ld1$modulus[1]
  ##ldet <- ld2 - ld1
  ldet <- ld1 - ld2
  
  ## rest of the calculation
  ## S2i <- solve(S2)
  S1i <- chol2inv(S1c)
  tr <- sum(diag(S1i %*% S2))
  m2mm1 <- mu2 - mu1
  qf <- as.numeric(t(m2mm1) %*% (S1i %*% m2mm1))

  ## return the correct combination of the parts
  r <- 0.5*(ldet + tr + qf - N)
  if(symm) return(0.5*(r + kl.norm(mu2, S2, mu1, S1, quiet, symm=FALSE)))
  else return(r)
}


## Ellik:
##
## the expected log likelihood under an MVN with parameters
## mu1 and S1 when supposing that the data follows an reference
## MVN distribution with paramerters mu2 and S2

Ellik.norm <- function(mu1, S1, mu2, S2, quiet=FALSE)
  {
    N <- length(mu1)
    
    ## check that the mean vectors have the same length
    if(length(mu2) != N) stop("must have length(mu1) == length(mu2)")
    
    ## check the covar matrices have same dims as the mean
    if(ncol(S1) != N || nrow(S1) != N)
      stop("must have nrow(S1) == ncol(S1) == length(mu1)")
    if(ncol(S2) != N || nrow(S2) != N)
      stop("must have nrow(S2) == ncol(S2) == length(mu2)")

    ## first calculate the differential entropy
    S2c <- chol(S2)
    ld2 <-  2*sum(log(diag(S2c)))
    de <- 0.5*(N*log((2*pi*exp(1))) + ld2)

    ## the calculate the KL-divergence
    S1c <- try(chol(S1), silent=TRUE)
    if(class(S1c) == "try-error") return(NA)
    ld1 <- 2*sum(log(diag(S1c)))
    S1i <- chol2inv(S1c)
    tr <- sum(diag(S1i %*% S2))
    m2mm1 <- mu2 - mu1
    qf <- as.numeric(t(m2mm1) %*% (S1i %*% m2mm1))
    ldet <- ld1 - ld2
    kl <- 0.5*(ldet + tr + qf - N)

    ## return the sum negated
    return(-(de+kl))
  }


## rmse.muS:
##
## calculate the Root Mean Squared Error between
## two MVN paramterizations (in terms of means
## and covariance matrices)

rmse.muS <- function(mu1, S1, mu2, S2)
  {
    resid.mu <- (mu1 - mu2)^2
    S1 <- as.vector(S1[upper.tri(S1, TRUE)])
    S2 <- as.vector(S2[upper.tri(S2, TRUE)])
    resid.S <- (S1 - S2)^2
    ##print(c(1, max(resid.mu), max(resid.S)))
    ##print(c(2, mean(resid.mu), mean(resid.S)))
    return(sqrt(mean(c(resid.mu, resid.S))))
  }
