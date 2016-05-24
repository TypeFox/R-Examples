# LatticeKrig  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2012
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2

LKDistComponents<- function(  x1, x2, delta, max.points = NULL, 
     mean.neighbor= 50,
     distance.type="Euclidean")
{
# find distances between x1 and x2 but only for those pairs where each 
# component is within a distance delta
#
   if( distance.type!="Euclidean"){
   	stop("Only Euclidean distance supported")
   	 	}
# figure out how large an array to allocate for distance matrix   
    n1 <- nrow(x1)
    n2 <- nrow(x2)
    dimension<- ncol(x1)
    if (is.null(max.points)) {
        Nmax <- max(c(n1,n2)) * mean.neighbor
    }
    else {
        Nmax <- max.points
    }
# note delta used here instead of delta^2  
   out <- .Fortran("LKdistComp", x1 = as.double(x1), n1 = as.integer(n1), 
                             x2 = as.double(x2), n2 = as.integer(n2),
                            dim = as.integer(dimension), 
                          delta = as.double(delta),
                            ind = as.integer(rep(0, Nmax * 2)),
                             rd = as.double(rep(-1, Nmax*dimension)),
                           Nmax = as.integer(Nmax),
                          iflag = as.integer(1),
                                PACKAGE = "LatticeKrig")
# NOTE rd is returned as a matrix 
    # negative iflag means one has run out of space
    if (out$iflag == -1) {
        cat("temp space set at", Nmax, fill = TRUE)
        stop("Ran out of space, increase the max.points")
    }
    # trim down to a sparse matrix object where the elements that have
    # nonzero values (there are out$Nmax of these)
    N <- out$Nmax
    # Note  distance information has nonzero indices in "spind" sparse matrix format
    # the component distances are organized as a matrix
    # each column of ra has the nonzero elements of a distance matrix 
    # corresponding to that column. 
    out<- list(ind = matrix(out$ind, nrow=Nmax, ncol=2         )[1:N,],
                ra = matrix(out$rd,  nrow=Nmax, ncol=dimension )[1:N,],
                da = c(n1, n2) )           
    return(out)
 }
