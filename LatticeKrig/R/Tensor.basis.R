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

Tensor.basis <- function(x1, centers, basis.delta,
                     max.points = NULL, 
                  mean.neighbor = 50,
                  BasisFunction = "WendlandFunction", 
                  distance.type = "Euclidean")
{    
    dimension <- ncol(x1)
    n1 <- nrow(x1)
    if (is.null(max.points)) {
        Nmax <- n1 * mean.neighbor
    }
    else {
        Nmax <- max.points
    }
    # evalute RBF basis functions at the x1 locations with
    # nodes given by centers. Returned is a list following the 
    # spind format but with a matrix in place of the vector of nonzero
    # matrix values 
    # This format is because a distance is required for each component
    # to evaluate the tensor product basis function
    out<- LKrigDistance( x1, centers,
                             delta =  basis.delta,
                        max.points = max.points,
                     mean.neighbor = mean.neighbor,
                     distance.type = distance.type,
                        components = TRUE)                    
    # evaluate distance  with tensor function on each coordinate
    # in this case ra is a matrix with each row being the componentwise
    # distances and with the maximum distance in any coordinate being
    # less than basis.delta.   
    out$ra <-  out$ra/basis.delta
    temp <- do.call( BasisFunction, list(d=out$ra[,1]))
    if (dimension > 1) {
        for (j in (2:dimension)) {
            temp <- temp * do.call( BasisFunction, list( d=out$ra[,j]))
        }
    } 
    out$ra<- temp  
    out <- spam(out[c("ind", "ra")], nrow=out$da[1], ncol= out$da[2] )
    return(out)
}
