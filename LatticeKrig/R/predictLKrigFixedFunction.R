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

predictLKrigFixedFunction <- function(object, xnew=NULL, Znew = NULL, drop.Z = FALSE){
     if( is.null(xnew)){
       xnew<- object$x
     }
     nt<- object$nt
     nZ<- object$nZ
# logical that indicates the spatial drift component (e.g. a low order polynomial)     
     ind.drift<- c( rep( TRUE, (nt-nZ) ), rep( FALSE, nZ)) 
#    cat( "predictLKrigFixedFunction:  ind.drift", ind.drift, fill=TRUE)
# predictions for fixed part of the model
# and can be with or without the additional covariates, Znew.
     distance.type<- object$LKinfo$distance.type
     T.matrix<- do.call(object$LKinfo$fixedFunction,
                            c(  list(x = xnew, Z = Znew,
                                  distance.type = distance.type),
                                  object$LKinfo$fixedFunctionArgs))
#     cat( "predictLKrigFixedFunction:  dim(T.matrix)", dim(T.matrix), fill=TRUE)                              
    if( !drop.Z){
      temp1<- T.matrix%*%object$d.coef}
    else{
      temp1<- T.matrix%*%object$d.coef[ind.drift, ]
    }
    return( temp1)
  }


