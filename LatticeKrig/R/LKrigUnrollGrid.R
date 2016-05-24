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

LKrigUnrollZGrid<- function( grid.list, ZGrid=NULL){
  if( is.null(ZGrid)){
    return(ZGrid)
  }
  if( is.list( ZGrid) ){
     if( any(grid.list[[1]] != ZGrid[[1]]) |any(grid.list[[2]] != ZGrid[[2]]) ){
         stop("grid list does not match grid for covariates")
       }  
# wipe out the x and y components of ZGrid
  ZGrid<- ZGrid$z
  }
# check dimensions
    Zdim<- dim( ZGrid)
      nx<- length( grid.list[[1]])
      ny<- length( grid.list[[2]])
      if( (Zdim[1] != nx) | (Zdim[2] != ny) ){
         stop( "Dimension of ZGrid does not match dimensions of location grid list.")
      }
# reshape as a matrix where rows index locations.
# Note that this works whether Zdim[3] exists or not! 
      return( matrix( c(ZGrid),  nrow= Zdim[1]*Zdim[2] ))
 }
  
