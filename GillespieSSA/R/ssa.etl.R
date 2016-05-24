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

`ssa.etl` <-
function(a = stop("missing propensity vector (a)"), 
        nu = stop("missing state-change matrix (nu)"),
       tau = stop("missing step size (tau)")) {
  M <- length(a)
  k <- rpois(M,(a*tau))
  # MPK: Strictly speaking it is not correct to call the realized state-change
  # vector nu_j here since, in leap methods, actually is not reaction specific.
  # In Pineda-Krch (JSS ms) it is refered to as the \bm{\tilde{nu}}
  return(list(tau=tau, nu_j=rowSums(matrix(rep(k,dim(nu)[1]),byrow=TRUE,ncol=M)*nu)))
}

