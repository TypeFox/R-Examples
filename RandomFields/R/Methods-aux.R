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


reflection <- function(data, orth, drop=FALSE)
  ##IMPORTANT NOTE! DO NOT CHANGE THE VARIABLE NAMES IN THIS SIGNATURE
  ## why ???
  ## since the variable data is pasted by its name
{
  d <- dim(data)
  return(do.call("[", c(list(data), rep(TRUE, orth-1), list(d[orth]:1),
                        rep(TRUE, length(d) - orth), drop=drop)))
}

AddUnits <- function(params) {
  ## see also empvario.R and fitgauss.R, if changed
  coords <- RFoptions()$general
  return(c(params, list(coordunits=coords$coordunits,
                        varunits=coords$varunits)))
}

compareGridBooleans <- function(grid, gridtmp) {
  if (!missing(grid) && length(grid)>0 && grid!=gridtmp)
    message(paste("you specified grid=", as.character(grid),
                  " but isGridded(data)=", as.character(gridtmp),
                  ";  grid is set to ", as.character(gridtmp), sep=""))
}

isSpObj <- function(x)
  (is(x, "SpatialGridDataFrame") || is(x, "SpatialPointsDataFrame")) &&
  !is(x, "RFsp")


sp2RF <- function(sp, param=list(n=1, vdim=1)) {
  class(sp) <- paste("RF", tolower(substr(class(sp), 1, 1)),
                       substring(class(sp), 2),  sep="")
  sp@.RFparams <- AddUnits(param)
  validObject(sp)
  return(sp)
}

convert2GridTopology <- function(grid){
  if (!is(grid, "GridTopology")) {
    if (is.null(dim(grid)))
      grid <- matrix(grid, ncol=1)
    stopifnot(nrow(grid)==3)
    grid <- sp::GridTopology(cellcentre.offset=grid[1,],
                             cellsize=grid[2,],
                             cells.dim=grid[3,])
  }
  return(grid)
}
     


## Generate Objects ########################################################

RFspatialGridDataFrame <- function(grid, data,
                                   proj4string = sp::CRS(as.character(NA)),
                                   RFparams=list(n=1, vdim=1)) {
  
  grid <- convert2GridTopology(grid)
  tmp <- sp::SpatialGridDataFrame(grid=grid,
                                  data = if (is.data.frame(data)) data else
                                  data.frame(data),
                                  proj4string=proj4string)
  return(sp2RF(tmp, RFparams))
#  tmp <- as(tmp, "RFspatialGridDataFrame")
#  tmp@.RFparams <- AddUnits(RFparams)
#  validObject(tmp)
#  return(tmp)
}

RFspatialPointsDataFrame <- function(coords, data, coords.nrs = numeric(0),
                                     proj4string = sp::CRS(as.character(NA)), 
                                     match.ID = TRUE, bbox = NULL,
                                     coordunits = NULL,
                                     varunits = NULL,
                                     RFparams=list(n=1, vdim=1)) {
  if (is.null(bbox)) {
    bbox <- t(apply(coords, 2, range))
    colnames(bbox) <- c("min", "max")    
  }
 
  tmp <- sp::SpatialPointsDataFrame(coords=coords,
                                    data=if (is.data.frame(data)) data else
                                    data.frame(data),
                                    coords.nrs=coords.nrs,
                                    proj4string=proj4string, 
                                    match.ID=match.ID, bbox=bbox)

  RFparams$n <- as.integer(RFparams$n)
  RFparams$vdim <- as.integer(RFparams$vdim)
  return(sp2RF(tmp, RFparams))

}

RFgridDataFrame <- function(data, grid,
                            RFparams=list()){
  grid <- convert2GridTopology(grid)
  data <- as.data.frame(data)
  return(new("RFgridDataFrame", data=data, grid=grid,
             .RFparams=AddUnits(RFparams)))
}

RFpointsDataFrame <- function(data=data.frame(NULL), coords=as.numeric(NULL),
                              RFparams=list()){
  data <- as.data.frame(data)
  if (is.null(dim(coords))) coords <- matrix(coords)
  return(new("RFpointsDataFrame", data=data, coords=coords,
             .RFparams=AddUnits(RFparams)))
}



brack <- function(x, i, j, ..., drop=FALSE) {
  dots = list(...)

  #Print(if (missing(x)) "no x" else x, if (missing(i)) "no i" else i,
  #      if (missing(j)) "no j" else j, dots)
  
  if (length(dots)>0) warning("dots are ignored")
  has.variance <- !is.null(x@.RFparams$has.variance) && x@.RFparams$has.variance
  if (missing(j)) {
    if (missing(i)) return(x)
    x@data <- x@data[i]#, drop=drop]
    n <- x@.RFparams$n
    v <- x@.RFparams$vdim
    if (!is.numeric(i)) {
      if (is.logical(i)) {
        i <- which(i)
      } else {
        stopifnot(all(i %in% colnames(x@data)))
        i <- match(i, colnames(x@data))
      }
    }
    if (length(unique(table(i%%v, rep(0, length(i))))) !=1 )
      stop(paste("for each variable selected, the same number of repetitions ",
                 "must be selected; you selected columns ",
                 paste(i, collapse=","), " but vdim=",v," and n=",n, sep=""))
    x@.RFparams$vdim <- v.new <- length(unique(i%%v))
    if (ret.has.var <- has.variance && any(i > n*v))
      x@.RFparams$has.variance <- ret.has.var
    x@.RFparams$n <- length(i) / v.new - ret.has.var
    
  } else {
    if(missing(i))  x@data <- x@data[,j]
    else x@data <- x@data[i,j]
  }
  return(x)
}


brack2 <- function(x, i, j, ..., value) {
  dots = list(...)
  if (length(dots)>0) warning("dots are ignored")
  if (missing(j)) 
    x@data[i] <- value
  else
    x@data[i,j] <- value
  return(x)
}



cbind_RFsp <- function(...) {  ##copied from sp package
  stop.ifnot.equal = function(a, b) {
    res = all.equal(a@grid, b@grid)
    if (!is.logical(res) || !res)
      stop("grid topology is not equal")
  }
  grds = list(...)
  ngrds = length(grds)
  if (ngrds < 1)
    stop("no arguments supplied")
  if (ngrds == 1)
    return(grds[[1]])
  ## verify matching topology:
  sapply(grds[2:ngrds], function(x) stop.ifnot.equal(x, grds[[1]]))
  gr = grds[[1]]
  gr@data = do.call(base::cbind, lapply(grds, function(x) x@data))
  ##for (i in 2:ngrds)
  ##	gr@data = cbind(gr@data, grds[[i]]@data)
  if (is(gr, "RFspatialGridDataFrame"))
    sp::proj4string(gr) = sp::CRS(sp::proj4string(grds[[1]]))
  gr
}

cbind_RFspPoints <- function(...) {  ##copied from sp package
  stop.ifnot.equal = function(a, b) {
    res = all.equal(a@coords, b@coords)
    if (!is.logical(res) || !res)
      stop("coords are not equal")
  }
  grds = list(...)
  ngrds = length(grds)
  if (ngrds < 1)
    stop("no arguments supplied")
  if (ngrds == 1)
    return(grds[[1]])
  ## verify matching topology:
  sapply(grds[2:ngrds], function(x) stop.ifnot.equal(x, grds[[1]]))
  gr = grds[[1]]
  gr@data = do.call(base::cbind, lapply(grds, function(x) x@data))
  ##for (i in 2:ngrds)
  ##	gr@data = cbind(gr@data, grds[[i]]@data)
  gr
}



extract.names <- function(names) {
  if (length(names) == 1) return(as.vector(names))
  nr <- strsplit(names[,1], ".")
  if (any(sapply(nr, length) != 2)) nr <- names[,1]
  else nr <- sapply(nr, function(x) x[1])

  nc <- strsplit(names[1,], ".")
  if (any(sapply(nc, length) != 2)) nc <- names[1,]
  else nc <- sapply(nc, function(x) x[1])

  return(list(nr, nc))
}



## Coerce Objects #########################################################
spatialGridObject2conventional <- function(obj, data.frame=FALSE) {
  timespacedim <- length(obj@grid@cells.dim)
  data <- as.matrix(obj@data)
  names <- colnames(data)
  
  has.variance <- !is.null(obj@.RFparams$has.variance) &&
    obj@.RFparams$has.variance
  dim(data) <- NULL
  vdimn <- c(obj@.RFparams$vdim, obj@.RFparams$n + has.variance)

  dim(data) <- c(obj@grid@cells.dim, vdimn)
#  Print(names)
  
  if (timespacedim > 1) data <- reflection(data, 2, drop=FALSE)
  ## re-ordering of 2nd space dimension since in sp objects, the 2nd dimension
  ## is in decreasing order


  if (data.frame) {
    dim(data) <- c(prod(obj@grid@cells.dim), prod(vdimn))
    colnames(data) <- names
    return(as.data.frame(data))
  }
  
  dim(names) <- vdimn
  vdim_close_together <- FALSE
  if (vdim_close_together) {
    perm <- c(timespacedim+1, 1:timespacedim, timespacedim+2) 
    data <- aperm(data, perm=perm)
    names <- aperm(names, perm[-1]) ### ?????
  }
  ## new order of dimensions: vdim, space-time-dims, n

  is.dim <- dim(data) != 1
  if (sum(is.dim) > 1) {
    dim(data) <- dim(data)[is.dim] # drop dimensions length 1
    l <- list()
    l[length(obj@grid@cells.dim) + (1:2)] <- extract.names(names)
    dimnames(data) <- l[is.dim]
  } else {
    dim(data) <- NULL
    #names(data) <- names
  }

  x <- rbind(obj@grid@cellcentre.offset,
             obj@grid@cellsize,
             obj@grid@cells.dim)

 # Print(obj, "TTTT", is(obj, "RFsp"))
  
  if (dimensions(obj)==1 ||
      !(ZF_GENERAL_COORD_NAME[2] %in% names(obj@grid@cellcentre.offset)))
    T <- NULL
  else {
    idxT1 <- which(ZF_GENERAL_COORD_NAME[2] ==names(obj@grid@cellcentre.offset))
    T <- x[,  idxT1]
    x <- x[, -idxT1, drop=FALSE]
  }

  .RFparams <- obj@.RFparams
  
  return(list(data=data, x=x, T=T, .RFparams=.RFparams, .names=names))
}

as.data.frame.RFpointsDataFrame <-
  as.data.frame.RFspatialPointsDataFrame <- function(x, ...) {
  #str(x); kkkk
  cbind(x@data, x@coords)
}
as.data.frame.RFgridDataFrame <-
  as.data.frame.RFspatialGridDataFrame <- function(x, ...) 
  spatialGridObject2conventional(x, TRUE)

setAs("RFspatialPointsDataFrame", "data.frame",
      function(from, to) from@data
)
setAs("RFspatialGridDataFrame", "data.frame",
      function(from, to) spatialGridObject2conventional(from, TRUE)$data
)



spatialPointsObject2conventional <- function(obj) {
  data <- as.matrix(obj@data)
  Enames <- names <- colnames(data)

  has.variance <-
    !is.null(obj@.RFparams$has.variance) && obj@.RFparams$has.variance
  dim(data) <- NULL
  vdimn <- c(obj@.RFparams$vdim, obj@.RFparams$n + has.variance)
  dim(data) <- c(nrow(obj@data), vdimn)  
  dim(Enames) <- vdimn
  Enames <- extract.names(Enames)
  vdim_close_together <- FALSE
  if (vdim_close_together) {
    perm <- c(2,1,3)
    data <- aperm(data, perm=perm)
    Enames <- aperm(Enames, perm[-1]) ### ?????
  }

  x <- obj@coords
  dimnames(x) <- NULL
  idxT1 <- which(ZF_GENERAL_COORD_NAME[2] == colnames(obj@coords))

#  Print(obj,ZF_GENERAL_COORD_NAME, colnames(obj@coords), length(idxT1),
 #       length(obj@.RFparams$T))
#  print(dimensions(obj))
#  print(!(ZF_GENERAL_COORD_NAME[2] %in% colnames(obj@coords)))
#  oooooo
  
  if (dimensions(obj)==1 || length(idxT1) + length(obj@.RFparams$T) == 0) {
    T <- NULL
    is.dim <- dim(data) != 1
    if (sum(is.dim) > 1) {    
      dim(data) <- dim(data)[is.dim] # drop dimensions length 1
      dimnames(data) <- c(list(NULL), Enames)[is.dim]
    } else {
      dim(data) <- NULL
      ##names(data) <- names
    }  
  } else {
    if (length(idxT1) == 0) idxT1 <- dimensions(obj)
    dimdata <- dim(data)
#    Print(dimdata)

    stopifnot(length(idxT1) == 1 || length(dimdata) != dimensions(obj))
    RFparams <- obj@.RFparams
    RFparams$n <- 1
    
 #   str(x)
    #Print(class(x), x, idxT1, x[, idxT1], coords=unique(x[, idxT1]),
   #       data=double(length(unique(x[,idxT1]))), RFparams)
    rpdf <- RFpointsDataFrame(coords=unique(x[, idxT1]),
                              data=double(length(unique(x[,idxT1]))),
                              RFparams=RFparams)
    T <- sp::points2grid(rpdf)
    
    if (obj@.RFparams$vdim==1) {
      dim(data) <- c(dimdata[1]/T@cells.dim, T@cells.dim, dimdata[-1:-2])
      dimnames(data) <- list(NULL,
                             paste("T", 1:T@cells.dim, sep=""),
                             Enames[[2]])
    } else {
      dim(data) <- c(dimdata[1], dimdata[2]/T@cells.dim,
                     T@cells.dim, dimdata[-1])
      dimnames(data) <- list(NULL,
                             paste("T", 1:T@cells.dim, sep=""),
                             Enames[[1]], Enames[[2]])
     
    }
    x <- x[1:(nrow(x)/T@cells.dim), -idxT1, drop=FALSE]
    T <- c(T@cellcentre.offset, T@cellsize, T@cells.dim)
  }

#  Print(data=data, x=x, T=T, .RFparams=obj@.RFparams)
  
  return(list(data=data, x=x, T=T, .RFparams=obj@.RFparams))
}


## convert 'RFsp' objects to conventional format of 'RFsimulate',
## i.e. data is an array and x a matrix of coordinates or gridtriple defs.

setGeneric(name = "RFspDataFrame2conventional", 
           function(obj, ...) standardGeneric("RFspDataFrame2conventional"))
setMethod("RFspDataFrame2conventional",
          signature=c("RFspatialGridDataFrame"),
          definition=spatialGridObject2conventional)
setMethod("RFspDataFrame2conventional", signature=c("RFgridDataFrame"),
          definition=spatialGridObject2conventional)
setMethod("RFspDataFrame2conventional",
          signature=c("RFspatialPointsDataFrame"),
          definition=spatialPointsObject2conventional)
setMethod("RFspDataFrame2conventional", signature=c("RFpointsDataFrame"),
          definition=spatialPointsObject2conventional)




prepare4RFspDataFrame <- function(model=NULL,
                                  info, x, y, z, T, grid=NULL, data, RFopt) {
  
  vdim <- info$vdim
  locinfo <- info$loc

#  Print(info)
  
  names <- GetDataNames(model=model, locinfo=locinfo,
                        coords=if (missing(x)) NULL else 
                        list(x=x,y=y, z=z, T=T, grid=grid))
   
  if (!is.null(names$varnames)) {
    if (vdim != length(names$varnames))
      stop(paste("you passed a formula for 'model' with left-hand side '",
                 paste(names$varnames, collapse=","),
                 "', but vdim of the model equals ", vdim, sep=""))
  }
  
  coordnames.incl.T <- names$coordnames
  
  ## coords or GridTopology 
  if (locinfo$grid) {
    coords <- NULL
    xgr <- cbind(locinfo$xgr, locinfo$T)
    colnames(xgr) <- coordnames.incl.T
    xgr[is.na(xgr)] <- 0
    gridTopology <- sp::GridTopology(xgr[1, ], xgr[2, ], xgr[3, ])
  } else { ## grid == FALSE
    gridTopology <- NULL
    
    # cbind of locations from x-matrix and T (if given)
    coords <- as.matrix(apply(t(locinfo$x), 2, rep,
                              times=(locinfo$totpts/locinfo$spatialtotpts)))
    if (locinfo$Zeit) {
      T <- locinfo$T
      coords <- cbind(coords, rep(seq(T[1], by=T[2], len=T[3]),
                                each=locinfo$spatialtotpts))
    }
    if (is.matrix(coords)) colnames(coords) <- coordnames.incl.T
  }

  if (RFopt$general$printlevel>=PL_IMPORTANT && RFopt$internal$warn_newstyle) {
    RFoptions(internal.warn_newstyle = FALSE)
    message("New output format of RFsimulate: S4 object of class 'RFsp';\n",
            "for a bare, but faster array format use 'RFoptions(spConform=FALSE)'.")
  }

  return(list(coords=coords, gridTopology=gridTopology, vdim=vdim, names=names))
}

### ist keine Methode im engeren Sinne. Habe ich aus Methods-RFsp.R
### rausgenommen, da bei jeglicher Aenderung in Methods-RFsp.R ich
### komplett neu installieren muss. Bei rf.R muss ich es nicht.
conventional2RFspDataFrame <-
  function(data, coords=NULL, gridTopology=NULL, n=1, vdim=1, T=NULL,
           vdim_close_together) {
  
  if (!xor(is.null(coords), is.null(gridTopology)))
    stop("one and only one of 'coords' and 'gridTopology' must be NULL")
  
  varnames <- attributes(data)$varnames
  ## may be NULL, if called from 'RFsimulate', the left hand side of model, if
  ## model is a formula, is passed to 'varnames'
  attributes(data)$varnames <- NULL
  
  ## grid case
  if (length(coords) == 0) {# war is.null(coords) -- erfasst coords=list() nicht
    grid <- convert2GridTopology(gridTopology) 
    timespacedim <- length(grid@cells.dim)

    ## naechste Zeile eingefuegt !! und (Martin 30.6.13) wieder
    ## auskommentiert. s. Bsp in 'RFsimulate'
    ## if (!is.null(dim(data)) && all(dim(data)[-1]==1)) data <- as.vector(data)
     
    if (is.null(dim(data))) {
      
      d <- c(grid@cells.dim,  if (vdim > 1) vdim,  if (n > 1) n)
      stopifnot(length(data) == prod(d))
      dim(data) <- d
      
      
    # stopifnot(1 == timespacedim + (n > 1) + (vdim > 1))               
      
    } else {
     
      if (length(dim(data)) != timespacedim + (n>1) + (vdim > 1)){          
        stop(paste(length(dim(data)),
                   "= length(dim(data)) != timespacedim + (n>1) + (vdim>1) =",
                   timespacedim, '+', (n>1), '+', (vdim > 1)))
      }
    }
        
    if (vdim>1 && vdim_close_together){
      ## new order of dimensions: space-time-dims, vdim, n
      perm <- c( 1+(1:timespacedim), 1, if (n>1) timespacedim+2 else NULL)
      data <- aperm(data, perm=perm)
    }
    if (timespacedim==1)
      call <- "RFgridDataFrame"
    else {
      ## 3/2015: unclear what to do if 1d space and time: also reflection??
      data <- reflection(data, 2, drop=FALSE)
      call <- "RFspatialGridDataFrame"
    }
  }
  
  ## coords case
  if (is.null(gridTopology)){
    if (vdim>1 && vdim_close_together){
      n.dims <- length(dim(data))
      perm <- c(2:(n.dims - (n>1)), 1, if (n>1) n.dims else NULL)
      data <- aperm(data, perm=perm)
    }
    if (is.null(dim(coords)) || ncol(coords)==1)
      call <- "RFpointsDataFrame"
    else call <- "RFspatialPointsDataFrame"
  }

  
  
  ## in both cases:
  dim(data) <- NULL
  data <- as.data.frame(matrix(data, ncol=n*vdim))
  
  if (is.null(varnames))
    varnames <- paste("variable", 1:vdim, sep="")
  if (length(varnames) == n*vdim)
    names(data) <- varnames
  else
    if (length(varnames) == vdim)
      names(data) <- paste(rep(varnames, times=n),
                           if (n>1) ".n", if (n>1) rep(1:n, each=vdim),sep="")
    else names(data) <- NULL
  
  ## Print(call, varnames, names(data), coords, is.null(coords))

  if (is.null(coords)){
#   Print(call, data=data, grid, RFparams=list(n=n, vdim=vdim, T=T))

    do.call(call, args=list(data=data, grid=grid,
                      RFparams=list(n=n, vdim=vdim, T=T)))
  } else {
    ##   Print(call, args=list(data=data, coords=coords,RFparams=list(n=n, vdim=vdim, T=T)))
    
     do.call(call, args=list(data=data, coords=coords,
                       RFparams=list(n=n, vdim=vdim, T=T)))
  }
}


