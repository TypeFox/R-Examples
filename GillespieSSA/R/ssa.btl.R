# Copyright 2007, 2008, 2010 Mario Pineda-Krch.
#
# This file is part of the R package GillespieSSA.
#
# GillespieSSA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# GillespieSSA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with GillespieSSA.  If not, see <http://www.gnu.org/licenses/>.

`ssa.btl` <-
function(x = stop("missing state vector (x)"),
         a = stop("missing propensity vector (a)"),
        nu = stop("missing state-change matrix (nu)"),
         f = stop("missing coarse-graining factor (f)")) { # See p.4 in Chatterjee et al. (2005)

  coercing <- FALSE

  # Calculate tau
  tau <- f/sum(a)   # Calculate the current tau
  if (tau>1) tau<-1 # Tau cannot be larger than unity!

  M <- length(a)    # Number of reaction channels
  tilde_x <- x    
  nu_j <- matrix(rep(0,length(x)))
 
  # Loop over all reaction channels having propensity fun>0 
  for (j in seq(M)[a>0]) {    
    if (any(nu[,j]<0)) { # do this if there are limiting reactions
      mask <- nu[,j]<0
      L <- min(floor(tilde_x[mask]/abs(nu[mask,j])))
      if (a[j]*tau>L) {
        p <- 1
        coercing <- TRUE
      }
      else p <- a[j]*tau/L  
      k <- rbinom(1,L,p)
    } 
    else { # do this if there are no limiting reactions
      k <- rpois(1,(a[j]*tau))
    }

    # Update tilde_x for the current reaction j
    tmp_nu_j <- matrix(rep(k,dim(nu)[1]), byrow=TRUE, ncol=1)*nu[,j]
    tilde_x <- tilde_x + tmp_nu_j

    # Record the current state change (it is returned by ssa.btl)
    nu_j <- nu_j + tmp_nu_j    
  } # for()

  # Throw a warning message if p was coerced to unity. Coercing implies too 
  # large steps-size due to a high coarse-graining factor (f)
  if(coercing) warning("coerced p to unity - consider lowering f")

  return(list(tau=tau, nu_j=nu_j))
}

