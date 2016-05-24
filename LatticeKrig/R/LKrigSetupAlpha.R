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

LKrigSetupAlpha <- function(object, ...){
  UseMethod("LKrigSetupAlpha")
}

LKrigSetupAlpha.default<- function( object, ...){
   alpha <- object$alpha
   nlevel <- object$nlevel
# some obvious defaults for alpha to simplify calls from LatticeKrig or LKrig.   
   if( is.na(alpha[1]) ) {
       if( nlevel==1){
          alpha<- 1.0
        }
       else{
# define alpha based on the power nu.       	
       	if( !is.null( object$nu)){
       		alpha<- exp( (1:nlevel)*object$nu )
       		alpha<- alpha/ sum( alpha)
       	}
       	else{
          alpha<- rep( NA, nlevel)
          }
       }
    }
    scalar.alpha <- !is.list(alpha) 
    if (scalar.alpha & (nlevel != 1) & (length(alpha) == 1)){
                stop( "Only one alpha specifed for multiple levels")}     
# coerce alpha to a list if it is passed as something else
    if (!is.list(alpha)) {
        alpha <- as.list(alpha)
    }
    return( alpha )
 }




