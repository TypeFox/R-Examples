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
                         
setGeneric(name="LKrigDistance", function( x1, x2, delta, ...  )
 { 
 	standardGeneric("LKrigDistance")
 	}
 )
    
setMethod("LKrigDistance", 
signature(x1 = "matrix", x2 = "matrix", delta = "numeric" ),
function( x1, x2, delta, 
         max.points = NULL, mean.neighbor = 50 ,  distance.type = "Euclidean",
         components = FALSE ){    
	if( !components){
		LKDist(x1, x2, delta,
		 max.points = max.points, mean.neighbor = mean.neighbor,
		 distance.type = distance.type)
		}
	else{
		LKDistComponents(x1, x2, delta,
	     max.points = max.points, mean.neighbor = mean.neighbor,
	     distance.type = distance.type)	
		}
	}                   
# end definition of method 
 )


setMethod("LKrigDistance", 
signature(x1 = "matrix", x2 = "gridList", delta = "numeric"),
function( x1, x2, delta, 
         max.points = NULL,    mean.neighbor = 50 ,  distance.type = "Euclidean",
         components = FALSE){    
	if( !components){
		LKDistGrid(x1, x2,  delta,
		 max.points = max.points, mean.neighbor = mean.neighbor,
		 distance.type = distance.type)	
		}
	else{
		LKDistGridComponents(x1, x2, delta,
	     max.points = max.points, mean.neighbor = mean.neighbor,
	     distance.type = distance.type)	
		}
	}                   
# end definition of method 
 )

     
 