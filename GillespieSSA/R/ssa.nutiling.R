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

ssa.nutiling <- function(a,nu,j) {
  M  <- dim(nu)[2]        # Number of reaction channels in nu-tile
  N  <- dim(nu)[1]        # Number of states in nu tile
  U  <- length(a)/M       # Number of tessallations of nu tile
  f  <- ceiling((j/M)-1)  # Frameshift factor
  jp <- j-f*M             # Relative reaction channel index
  nu_jp <- nu[,jp] 
  nu_j <- c(rep(0,f*N),   # Leading zeros 
            nu_jp,        # Relative state-change matrix
            rep(0,(U*N-(f*N+N)))) # Lagging zeros
  return(nu_j)
}
