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


## randmvn:
##
## returns N random (x) samples from a randomly generated normal
## distribution via (mu) and (S), also returned.  mu is a
## standard normal vector of d means, and S is an inverse-Wishart
## renerated (dxd) covariance matrix with d+2 degrees of freedom
## and an identity centering matrix (mean)

'randmvn' <-
function(N, d, method=c("normwish", "parsimonious"),
         mup=list(mu=0, s2=1), s2p=list(a=0.5, b=1),
         pnz=0.1, nu=Inf)
  {
    ## check N
    if(length(N) != 1 || N < 0)
      stop("N should be a nonnegative integer")

    ## check d
    if(length(d) != 1 || d <= 0)
      stop("d should be a positive integer")

    ## check mup
    if(length(mup$mu) != 1)
      stop("mup$b should be a scalar")
    if(length(mup$s2) != 1 || mup$s2 <= 0)
      stop("mup$s2 should be a positive scalar")
    
    ## check s2p
    if(length(s2p$b) != 1 || s2p$b < 0)
      stop("s2$b should be a nonnegative scalar")
    if(length(s2p$a) != 1 || s2p$a <= 0)
      stop("s2p$a should be a positive scalar")

    ## check method
    method <- match.arg(method)

    ## check nu
    if(length(nu) != 1 || d < 1)
      stop("nu should be a scalar >= 1")

    ## using normal and inverse wishart method
    if(method == "normwish") {
    
      ## generate a coavariance matrix and mean vector
      s2 <- 1.0##/rgamma(1, shape=s2p$a, scale=1/s2p$b)
      S <- solve(rwish(d+2, s2*diag(d)))
      mu <- rnorm(d, mup$mu, mup$s2)
      
    } else if(method == "parsimonious") {
      ## using the beta set to zero method

      ## initialize empty mu and S
      mu <- rep(NA, d)
      S <- matrix(NA, nrow=d, ncol=d)

      ## first vector requires no regressions
      mu[1] <- rnorm(1, mup$mu, sqrt(mup$s2))
      S[1,1] <- 1.0/rgamma(1, shape=s2p$a, scale=1/s2p$b)

      ## use monomvn style regressions for the rest of the columns
      for(i in 2:d) {
        icept <- rnorm(1, mup$mu, mup$s2)
        alpha <-  s2p$a+i-1
        s2 <- 1.0/rgamma(1, shape=alpha, scale=1/s2p$b)
        beta <- rnorm(i-1, sd=sqrt(s2))
        nz <- rbinom(1, i-2, 1-pnz)+1 ## sample(1:(i-1), 1)
        beta[sample(1:(i-1), nz)] <- 0
        mu[i] <- icept + t(beta) %*% mu[1:(i-1)]
        S[i,1:(i-1)] <- t(beta) %*% S[1:(i-1),1:(i-1)]
        S[i,i] <- s2 + S[i,1:(i-1)] %*% beta
        S[1:(i-1),i] <- S[i,1:(i-1)]
      }
    }

    ## don't draw if N=0
    if(N == 0) return(list(mu=mu, S=S))
      
    ## draw N samples from the MVN or Student-t
    if(is.infinite(nu)) x <- rmvnorm(N, mu, S)#, method="chol")
    else {
      x <- rmvt(N, S, nu)
      x <- t(t(x) + mu)
    }
    return(list(x=x, mu=mu, S=S))
  }
