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

LKrigSetupAwght <- function(object,...){
  UseMethod("LKrigSetupAwght")
}

LKrigSetupAwght.default<- function( object,...){
# object == LKinfo  
   a.wght<- object$a.wght
  nlevel<- object$nlevel
  if (!is.list(a.wght)) {
        # some checks on a.wght
        # coerce a.wght to list if it is passed as something
        # else (most likely a vector)
        if (nlevel == 1) {
            a.wght <- list(a.wght)
        }
        else {
            # repeat a.wght to fill out for all levels.
            if (length(a.wght) == 1) {
                a.wght <- rep(a.wght, nlevel)
            }
            a.wght <- as.list(c(a.wght))
        }
    }
   # check length of a.wght list
    if (length(a.wght) != nlevel) {
        stop("length of a.wght list differs than of nlevel")
    }
# default is to use the usual cholesky based normalization    
# see LKRectangle for an example of something different.
     attr( a.wght, "fastNormalize") <- FALSE
     return(a.wght)
 }
