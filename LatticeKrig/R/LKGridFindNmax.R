
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

LKGridFindNmax<- function(n1, max.points, mean.neighbor, delta, gridList){    
# figure out how large an array to allocate for distance matrix   
#    n1 <- nrow(x1)
    info<- summary( gridList)
    deltaScaled<- delta/info$dx
# estimate upper bound for the max.points
    if( !is.null( mean.neighbor) ){
	    max.points <- mean.neighbor*n1 
	 }   
    if (is.null(max.points)) {
        max.points <- n1 * ceiling(prod(deltaScaled*2 + 1 ) )
     }
    return( max.points )
} 
