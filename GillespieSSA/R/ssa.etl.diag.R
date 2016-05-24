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

`ssa.etl.diag` <- function(a, nu_tile, tau) {
  MU <- length(a)         # Toto nr of reaction channels
  k  <- rpois(MU,(a*tau)) # Nr of firings per channel
  M  <- dim(nu_tile)[2]   # Nr of reaction channel per patch (nu_tile)
  U  <- MU/M              # Nr of tilings
  nu_j <- NULL
  for(f in (seq(U)-1))
    nu_j <- c(nu_j, rowSums(matrix(rep(k[1:M+f*M],dim(nu_tile)[1]),byrow=TRUE,ncol=M)*nu_tile))
  return(list(tau=tau, nu_j=nu_j))
}

