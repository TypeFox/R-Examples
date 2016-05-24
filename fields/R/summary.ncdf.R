# fields  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2016
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

"summary.ncdf" <- function (object, ...) 
{
     tempList<- NULL
     varNames<- NULL
#
    cat("DIMENSIONS", fill=TRUE)
    for (i in names(object$dim) ) {
        vname = i
        ndims = length(object$dim[[i]]$vals)
         cat(vname, " has size", ndims, fill=TRUE)
    }
    cat(fill=TRUE) 
    cat("VARIABLES", fill=TRUE)
    for (i in 1:object$nvars) {
        vname = object$var[[i]]$name
        ndims = object$var[[i]]$ndims
        dimstring = paste(vname, "( variable ",i , 
            ") has shape")
        dimTemp<- NULL
        for (j in 1:ndims) {
            dimTemp<- c( dimTemp, object$var[[i]]$dim[[j]]$len)
        }
        temp<- ( dimTemp)
        varNames<- c(varNames, vname)
        tempList<- c( tempList, list(dimTemp))
        if( is.null(dimTemp) ){
           dimTemp<- NA}
        cat( i,":",  vname, "has size ", dimTemp, sep=" ", fill = TRUE)
    }
     names(tempList) <- varNames
     invisible( tempList)
}
