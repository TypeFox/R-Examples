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




setAs("SpatialPointsDataFrame", "RFspatialPointsDataFrame",
      function(from, to) sp2RF(from))
setAs("SpatialGridDataFrame", "RFspatialGridDataFrame",
      function(from, to) sp2RF(from))

setAs("RFspatialGridDataFrame", "RFspatialPointsDataFrame",
      function(from, to) RFspatialPointsDataFrame(coords=coordinates(from),
                                                  data=from@data,
                                                  RFparams=from@.RFparams))
setAs("RFspatialPointsDataFrame", "RFspatialGridDataFrame",
      function(from, to) RFspatialGridDataFrame(grid=sp::points2grid(from),
                                                data=from@data,
                                                RFparams=from@.RFparams))

setAs("RFpointsDataFrame", "RFgridDataFrame",
      function(from, to) {
        tmp.coord <- sp::SpatialPoints(cbind(from@coords, 1:length(from@coords)))
        tmp.grid <- sp::points2grid(tmp.coord)
        grid <- convert2GridTopology(tmp.grid[,1, drop=FALSE])
        RFgridDataFrame(data=from@data, grid=grid, RFparams=from@.RFparams)
      })

grid2pts1D <- function(from) {
        coords <- GridTopology2gridVectors(from@grid)[[1]]
        RFpointsDataFrame(data=from@data, coords=coords,
                          RFparams=from@.RFparams)
      }
setAs("RFgridDataFrame", "RFpointsDataFrame",
      function(from, to) grid2pts1D(from))



## methods 'as.matrix' and 'cbind' as in the 'sp' package

as.matrix.RFgridDataFrame <- 
 as.matrix.RFspatialGridDataFrame <- function(x, ...) {
  z <- as.array.RFgridDataFrame(x, ...) 
  if (length(dim(z)) > 2)
    stop("the data set cannot be turned into a matrix. Try 'as.array'")
  z
}

as.array.RFgridDataFrame <-
  as.array.RFspatialGridDataFrame <- function(x, ...) {
    z <- as.matrix(x@data)    
    dim <- c(x@grid@cells.dim, x@.RFparams$vdim, x@.RFparams$n)
    dim <- dim[dim > 1]
    if (length(dim) == 1) return(as.vector(z))
    dim(z) <- dim
    if ((l <- length(x@grid@cells.dim)) > 1) {
      idx <- as.list(rep(TRUE, length(dim)))
      idx[[2]] <- (dim[2]) : 1
      #Print(idx, c(list(z), drop=FALSE, idx))
      do.call("[", c(list(z), drop=FALSE, idx))
    } else z
  }


as.vector.RFgridDataFrame <-
  as.vector.RFspatialGridDataFrame <- function(x, ...) {
    as.vector(as.array.RFgridDataFrame(x))
  }

as.matrix.RFpointsDataFrame <-
  as.matrix.RFspatialPointsDataFrame <- function(x, ...) {
  z <- as.array.RFspatialPointsDataFrame(x, ...) 
  if (length(dim(z)) > 2)
    stop("the data set cannot be turned into a matrix. Try 'as.array'")
  z
}
as.array.RFpointsDataFrame <-
  as.array.RFspatialPointsDataFrame <- function(x, ...) {
  z <- as.matrix(x@data)
  dim <- c(x@.RFparams$vdim, x@.RFparams$n)
  if (length(x@.RFparams$T) == 3) dim <- c(x@.RFparams$T[3], dim)
  dim <- c(length(z) / prod(dim), dim)
  dim <- dim[dim > 1]
  if (length(dim) == 1) return(as.vector(z))
  dim(z) <- dim
  z
}
as.vector.RFpointsDataFrame <-
  as.vector.RFspatialPointsDataFrame <- function(x, ...) {
  as.vector(as.array.RFspatialPointsDataFrame(x))
}



#setAs("RFgridDataFrame", "matrix", gridtomatrix)
#setAs("RFspatialGridDataFrame", "matrix", gridtomatrix)


cbind.RFspatialGridDataFrame <- function(...)
  cbind_RFsp(...)
cbind.RFgridDataFrame <- function(...)
  cbind_RFsp(...)

cbind.RFspatialPointsDataFrame <- function(...)
  cbind_RFspPoints(...)
cbind.RFpointsDataFrame <- function(...)
  cbind_RFspPoints(...)

range.RFspatialGridDataFrame <- function(..., na.rm = FALSE, finite = FALSE)
  base::range(as.vector(...), na.rm = na.rm, finite = finite)
range.RFgridDataFrame <- function(..., na.rm = FALSE, finite = FALSE)
  base::range(as.vector(...), na.rm = na.rm, finite = finite)
range.RFspatialPointsDataFrame <- function(..., na.rm = FALSE, finite = FALSE)
  base::range(as.vector(...), na.rm = na.rm, finite = finite)
range.RFpointsDataFrame <- function(..., na.rm = FALSE, finite = FALSE)
  base::range(as.vector(...), na.rm = na.rm, finite = finite)


hist.RFspatialGridDataFrame <- function(x, ...)
  graphics::hist(as.vector(x), ...)
hist.RFgridDataFrame <- function(x, ...)
  graphics::hist(as.vector(x), ...)
hist.RFspatialPointsDataFrame <- function(x, ...)
  graphics::hist(as.vector(x), ...)
hist.RFpointsDataFrame <- function(x, ...)
  graphics::hist(as.vector(x), ...)


rfspDataFrame2dataArray_generic <- function(obj) {
  has.variance <- !is.null(obj@.RFparams$has.variance) && obj@.RFparams$has.variance
  arr <- unlist(obj@data)
  names(arr) <- NULL
  grid.length.v <- obj@grid@cells.dim
  dim(arr) <- c(grid.length.v,         # spatial and time dim
                obj@.RFparams$vdim,    # vdim
                obj@.RFparams$n + has.variance)       # repetitions
  timespacedim <- length(grid.length.v)
  if (timespacedim > 1)
    arr <- reflection(arr, 2, drop=FALSE)
    ## re-ordering of 2nd space dimension since in sp objects, the 2nd dimension
    ## is in decreasing order
  return(arr)
}      

setGeneric(name = "RFspDataFrame2dataArray", 
           def = function(obj) standardGeneric("RFspDataFrame2dataArray"))
setMethod("RFspDataFrame2dataArray", signature=c("RFspatialGridDataFrame"),
          definition=rfspDataFrame2dataArray_generic)
          ## used in plot method
setMethod("RFspDataFrame2dataArray", signature=c("RFgridDataFrame"),
          definition=rfspDataFrame2dataArray_generic)



##  convert GridTopology or 3-row matrix to list of vectors
setGeneric(name = "GridTopology2gridVectors", 
           def = function(grid) standardGeneric("GridTopology2gridVectors"))
setMethod("GridTopology2gridVectors",
          signature=c("GridTopology"),
          function(grid) {
            len <- length(grid@cells.dim)
            x <- list()
            for (i in 1:len)
              x[[i]] <- seq(from   = grid@cellcentre.offset[i],
                            to     = grid@cellcentre.offset[i] +
                                     ((grid@cells.dim[i]-1) * grid@cellsize[i]),
                            length = grid@cells.dim[i])
            return(x)
          })
setMethod("GridTopology2gridVectors",
          signature=c("matrix"),
          function(grid) {
            stopifnot(nrow(grid)==3)
            x <- list()
            for (i in 1:ncol(grid))
              x[[i]] <- seq(from  = grid[1,i],
                            to    = grid[1,i] + (grid[3,i]-1) * grid[2,i],
                            length= grid[3,i])
            return(x)
          })
##  convert GridTopology to 3-row matrix
as.matrix.GridTopology <- function(x, ...)
  return(rbind(x@cellcentre.offset, x@cellsize, x@cells.dim))



## extract or calculate! coordinates; 
                                                
setMethod(f = "coordinates", signature="RFpointsDataFrame",
          definition=function(obj) return(obj@coords) )
setMethod(f = "coordinates", signature="RFgridDataFrame",
          definition=function(obj) coordinates(as(obj, "RFpointsDataFrame")) )


extract.kriging.variance <- function(obj){
  stopifnot(is(obj, "RFsp"))
  has.variance <-
    !is.null(obj@.RFparams$has.variance) && obj@.RFparams$has.variance
  if (!has.variance) return(NULL)
  return(obj[(obj@.RFparams$n*obj@.RFparams$vdim + 1):(ncol(obj@data))])
}

  
setGeneric(name = "variance", 
           def = function(obj) standardGeneric("variance"))
setMethod(f = "variance", signature="RFspatialGridDataFrame",
          definition = extract.kriging.variance)
setMethod(f = "variance", signature="RFspatialPointsDataFrame",
          definition = extract.kriging.variance)
setMethod(f = "variance", signature="RFgridDataFrame", 
          definition = extract.kriging.variance)
setMethod(f = "variance", signature="RFpointsDataFrame",
          definition = extract.kriging.variance)



## conventional 'RFsimulate' output to 'RFsp' class

summary.RFpointsDataFrame <- function(object, digits = 6, ...) {
  df = data.frame(coordinates=signif(coordinates(object), digits), object@data)
  row.names(df) = row.names(object@data)
  class(df) <- "summary.RFpointsDataFrame"
  df
}
print.summary.RFpointsDataFrame <- function(x, ...) {
  if (is.data.frame(x)) print.data.frame(x, ...)
  else str(x, give.attr=FALSE, give.head=FALSE, no.list=TRUE, ...) #
}

print.RFpointsDataFrame <- function(x, ...) {
  sx <- summary.RFpointsDataFrame(x, ...)
  cat("Object of class RFpointsDataFrame\n")
  print.summary.RFpointsDataFrame(sx) #
  invisible(sx)
}

setMethod(f="show", signature="RFpointsDataFrame",
          definition=function(object) print.RFpointsDataFrame(object))

summary.RFgridDataFrame <- function(object, ...) {
  if (ncol(object@data) > 1) summary(object@data) else summary(object@data[[1]])
} 
print.RFgridDataFrame <- function(x, ...) {
  cat("Object of class RFgridDataFrame\n")
  cat("Grid topology:\n")
  print.data.frame(data.frame(cellcentre.offset=x@grid@cellcentre.offset, #
                   cellsize=x@grid@cellsize,
                   cells.dim=x@grid@cells.dim, row.names=""))
  cat("Points:\n")
  print.data.frame(data.frame(coordinates=coordinates(x), x@data)) #
  cat("Data summary:\n")
  summary.RFgridDataFrame(x, ...)
  invisible(x)
}

setMethod(f="show", signature="RFgridDataFrame", 
	  definition=function(object) print.RFgridDataFrame(object))

summary.RFspatialPointsDataFrame <- function(object, ...) summary(object@data)
print.RFspatialPointsDataFrame <- function(x, ...){
  if (!hasArg("silent") || !list(...)$silent) {
    cat("Object of class 'RFspatialPointsDataFrame'\n")
    str(x@data, no.list=TRUE, give.attr=FALSE)#
  }
  invisible(x@data)
}
setMethod(f="show", signature="RFspatialPointsDataFrame", 
	  definition=function(object) print.RFspatialPointsDataFrame(object))


summary.RFspatialGridDataFrame <- function(object, ...) summary(object@data)
print.RFspatialGridDataFrame <- function(x,...) {
 if (!hasArg("silent") || !list(...)$silent) {
    cat("Object of class 'RFspatialGridDataFrame'\n")
    utils::str(x@data, no.list=TRUE, give.attr=FALSE)#
  }
  invisible(x@data)
}
setMethod(f="show", signature="RFspatialGridDataFrame", 
	  definition=function(object) print.RFspatialGridDataFrame(object))


## extend methods for 'Spatial' to class 'RFsp'
if (!isGeneric("isGridded"))
  setGeneric(name = "isGridded", 
             def = function(obj) standardGeneric("isGridded") )

#setMethod(f="isGridded", signature="RFsp", 
#	  definition=function(obj)
#          (is(obj, "RFgridDataFrame") || is(obj, "RFspatialGridDataFrame")))
setMethod(f="isGridded", "RFgridDataFrame", function(obj) TRUE)
setMethod(f="isGridded", "RFpointsDataFrame", function(obj) FALSE)
setMethod(f="isGridded", "RFspatialGridDataFrame", function(obj) TRUE)
setMethod(f="isGridded", "RFspatialPointsDataFrame", function(obj) FALSE)
       
setMethod(f="dimensions", signature="RFspatialGridDataFrame", 
	  function(obj) (getMethod("dimensions", "Spatial")@.Data)(obj) )
setMethod(f="dimensions", signature="RFspatialPointsDataFrame", 
	  function(obj) (getMethod("dimensions", "Spatial")@.Data)(obj) )
setMethod(f="dimensions", signature="RFgridDataFrame", 
	  definition=function(obj) return(1))
setMethod(f="dimensions", signature="RFpointsDataFrame", 
	  definition=function(obj) return(1))

