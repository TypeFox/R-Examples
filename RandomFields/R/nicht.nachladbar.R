
## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2015 Martin Schlather
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 3
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  



brack2<- function(x, i, j, ..., value) {
            dots = list(...)
            if (length(dots)>0) warning("dots are ignored")
            if (missing(j)) 
              x@data[i] <- value
            else
              x@data[i,j] <- value
            return(x)
          }

setMethod("[", signature=c("RFgridDataFrame"), def=brack)
setMethod("[", signature=c("RFpointsDataFrame"), def=brack)
setMethod("[", signature=c("RFspatialGridDataFrame"), def=brack)
setMethod("[", signature=c("RFspatialPointsDataFrame"), def=brack)

setMethod("[<-", signature=c("RFgridDataFrame", "ANY", "ANY"),
          def=brack2)
setMethod("[<-", signature=c("RFpointsDataFrame", "ANY", "ANY"),
          def=brack2)
setMethod("[<-", signature=c("RFspatialGridDataFrame", "ANY", "ANY"),
          def=brack2)
setMethod("[<-", signature=c("RFspatialPointsDataFrame", "ANY", "ANY"),
          def=brack2)





trafo_pointsdata <- function(x, dim) {
  if (isgrid <- is(x, "RFgridDataFrame")) {
 #   Print(x)
    x <- grid2pts1D(x) ## x <- as(x, "RFpointsDataFrame") # funktionierte nicht
    ##    Print(x)
  } else if ((is(x, "matrix") || is(x, "data.frame")) && !missing(dim)) {
    dc <- data.columns(x, xdim = dim, force=TRUE)
    x <- list(coords=x[, dc$x, drop=FALSE], data=x[, dc$data, drop=FALSE])
  } else {
    if (!is(x, "RFpointsDataFrame"))
      stop("method only for objects of class 'RFpointsDataFrame' and 'RFgridDataFrame'")
  }
             
  dummy <- dimnames(x@coords)[[2]][1]
  lab <- xylabs(if (is.null(dummy)) "" else dummy, "",
                units=x@.RFparams$coordunits)
  labdata <- names(x@data)
  colname <- colnames(x@data)
  if (isgrid) {
    return(list(coords=as.vector(x@coords),
                data=as.matrix(x@data),
                RFparams=x@.RFparams,
                lab=lab, labdata=labdata, colnames=colname))
  } else {
    ord <- order(x@coords)
    return(list(coords=x@coords[ord, ],
                data=as.matrix(x@data)[ord, , drop=FALSE],
                RFparams=x@.RFparams,
                lab=lab, labdata=labdata, colnames=colname))
  } 
}
