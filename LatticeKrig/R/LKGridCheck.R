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

# A utility function to setup up some variables and avoid 
# dupicating code between radial and component functions 
 LKGridCheck<- function( distance.type, x1, gridList){
 	 if( distance.type!="Euclidean"){
     	stop("Only Euclidean distance supported")
   	 	}	
     if( !is.matrix( x1)){
     	stop( "x1 must be a matrix")
     }	  	 	
   	 if( class(gridList)!="gridList"){
   	 	stop("gridList must have class gridList")
   	 }
   	info<- summary( gridList)
   	if( info$dim != ncol(x1)){
   		stop("number of components in gridList and columns of x1
   		are not equal.")
   	}
# check that all the grid spacings are the same    
     check<- max( abs(info$dx - mean( info$dx))/ mean( info$dx) )
    if( check> 1e-7 )  {
    	stop(paste( check, " gridList spacings must all be the same.") )
    }
 }

