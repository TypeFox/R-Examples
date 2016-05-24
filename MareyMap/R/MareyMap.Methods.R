# Copyright 2006 Laboratoire de Biologie et de Biometrie Appliquée 
# (UMR 5558);CNRS; Univ. Lyon 1, 43 bd 11 nov, 69622, 
# Villeurbanne Cedex, France.
#
# This file is part of MareyMap.
#
# MareyMap is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# MareyMap is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MareyMap; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA



#--- constructors
setMethod("MareyMap", signature(data_table = "missing", column_names = "missing", set_name = "missing", map_name = "missing"), function() {
  new("MareyMap")
})


setMethod("MareyMap", signature(data_table = "data.frame"), function(data_table, column_names = list(name = "name", phys = "phys", gen = "gen", vld = "vld"), set_name = "", map_name = "") {
  object <- new("MareyMap")
  object@setName <- set_name
  object@mapName <- map_name
  object@markerNames <- as.vector(data_table[[column_names$name]])
  object@physicalPositions <- as.vector(data_table[[column_names$phys]])
  object@geneticDistances <- as.vector(data_table[[column_names$gen]])
  object@markerValidity <- as.vector(data_table[[column_names$vld]])
  object
})

  
#--- accessors
setMethod("setName", "MareyMap", function(object) object@setName)
setMethod("mapName", "MareyMap", function(object) object@mapName)
setMethod("markerNames", "MareyMap", function(object) object@markerNames)
setMethod("physicalPositions", "MareyMap", function(object) object@physicalPositions)
setMethod("geneticDistances", "MareyMap", function(object) object@geneticDistances)
setMethod("interpolations", "MareyMap", function(object) object@interpolations)
setMethod("markerValidity", "MareyMap", function(object) object@markerValidity)    


#--- replace methods
setReplaceMethod("physicalPositions", "MareyMap", function(object, value) {
  object@physicalPositions <- value
  object
})


setReplaceMethod("geneticDistances", "MareyMap", function(object, value) {
  object@geneticDistances <- value
  object
})


setReplaceMethod("markerNames", "MareyMap", function(object, value) {
  object@markerNames <- value
  object
})


setReplaceMethod("setName", "MareyMap", function(object, value) {
  object@setName <- value
  object
})


setReplaceMethod("mapName", "MareyMap", function(object, value) {
  object@mapName <- value
  object
})


setReplaceMethod("interpolations", "MareyMap", function(object, value) {
  object@interpolations <- value
  object
})


setReplaceMethod("markerValidity", "MareyMap", function(object, value) {
  object@markerValidity <- value
  object
})


#--- interpolationNames
setGeneric("interpolationNames", function(object, set_name, map_name) standardGeneric("interpolationNames"))


#--- interpolation
setMethod("interpolation", c("MareyMap", "character"), function(object, inter_name) {
  inter <- object@interpolations
  inter[match(inter_name, lapply(inter, function(x){name(x)}))][[1]]
})


#--- interpolation<-
setMethod("interpolation<-", c("MareyMap", "character", "Interpolation"), function(object, inter_name, value) {
  object@interpolations[match(inter_name, lapply(object@interpolations, function(x) {name(x)}))][[1]] <- value
  object
})


setMethod("+", c(e1 = "MareyMap", e2 = "Interpolation"), function(e1, e2) {
  e2 <- interpolate(e2, e1)
  e1@interpolations[[name(e2)]] <- e2
  e1
})


#---
setMethod("validPositions", "MareyMap", function(object) {
  physicalPositions(object)[which(markerValidity(object))]
})


#--- plot
setMethod("plot", signature(x = "MareyMap", y = "missing"), function(x , y,...) {
  cl <- match.call()
  cl$object <- cl$x
  cl$x <- NULL
  
	#split the screen and draw both                
  lo <- layout(matrix(c(1, 2), 2), heights = c(2, 1))
  layout.show(lo)
  
  if(is.null(cl$model) || cl$model) {
    cl[[1]] <- as.name("plotModels")
    par(mar = c(4, 5, 1, 1) + 0.1)
    if(is.null(cl$markers) || cl$markers) {
			# don't draw axis if the map is to be drawn as well
      cl$xaxt <- "n"
      cl$yaxt <- "n"
      cl$ylab <- ""        
      cl$xlab <- ""
    }
    eval(cl)
  }
  
  if(is.null(cl$markers) || cl$markers) {
    par(new = TRUE)
    cl[[1]] <- as.name("plotMarkers")
    cl$xaxt <- NULL
    cl$yaxt <- NULL
    cl$xlab <- NULL
    cl$ylab <- NULL
    cl$ylim <- NULL
    eval(cl)
	  cl[[1]] <- as.name("plotRates")
	  eval(cl)
  }
})

  
#--- plotMarkers
setMethod("plotMarkers", "MareyMap", function(object, ...) {
  par(mar = c(4, 5, 1, 1) + 0.1)
  cl <- match.call()
  if(is.null(cl$xlab))
    cl$xlab <- "Physical positions (Mb)" 
  if(is.null(cl$ylab))
    cl$ylab <- "Genetic distances (cM)"
  plot(object@physicalPositions / 1000000, object@geneticDistances, xlab = cl$xlab, ylab = cl$ylab, ...)
  if(length(which(!markerValidity(object))) != 0)
    points(physicalPositions(object)[which(!markerValidity(object))] / 1000000, geneticDistances(object)[which(!markerValidity(object))], col = "red", pch = 3, type = "p")
})


#--- plotModels
setMethod("plotModels", "MareyMap", function(object, ...) {
  par(mar = c(4, 5, 1, 1) + 0.1)
  cl <- match.call()
  if(is.null(cl$ylim))
    cl$ylim <- c(0, max(geneticDistances(object), na.rm = TRUE))
  if(is.null(cl$xlim))
    cl$xlim <- c(0, max(physicalPositions(object) / 1000000, na.rm = TRUE))
  if(is.null(cl$xlab))
    cl$xlab <- "Physical positions (Mb)"
  if(is.null(cl$ylab))    
    cl$ylab <- "Genetic distances (cM)"
  plot(-1, -1, xlab = cl$xlab, ylab = cl$ylab, col = "white", xaxt = "n", yaxt = "n")
  cl$xaxt <- "n"
  cl$yaxt <- "n"
  cl$xlab <- ""
  cl$ylab <- ""
  cl$col <- NULL
  cl$x <- NULL
  cl$y <- NULL
  cl[[1]] <- as.name("plotModel")
  lapply(interpolations(object), function(itr) {
    if(visible(itr)) {
      par(new = TRUE)
      cl$object <- itr
      eval(cl)
    }
  })
  return()
})


#--- plotRates
setMethod("plotRates", "MareyMap", function(object, ...) {
  cl <- match.call()
  par(mar = c(4, 5, 1, 1) + 0.1)
  if(length(object@interpolations) == 0) {
    par(col = "white", col.lab = "white", col.axis = "white", tck = 0)
    plot(0, 0)
  } else {
    if(is.null(cl$xlim))
      cl$xlim <- c(0, max(physicalPositions(object) / 1000000, na.rm = T))
    if(is.null(cl$ylim))
      cl$ylim <- c(0, max(unlist(lapply(object@interpolations, function(itr) {max(rates(itr), na.rm = TRUE)}))) + 1)
    if(is.null(cl$xlab))
      cl$xlab <- "Physical positions (Mb)"
    if(is.null(cl$ylab))
      cl$ylab <- "Recomb. rates \n (cM/Mb)"
		## setting up axis
    plot(0, 0, xlab = cl$xlab, ylab = cl$ylab, xlim = cl$xlim, ylim = cl$ylim, col = "white")
    
    cl$xaxt <- "n"
    cl$yaxt <- "n"
    cl$ylab <- ""        
    cl$xlab <- ""
    cl[[1]] <- as.name("plotRate")
    lapply(interpolations(object), function(itr) {
      if(visible(itr)) {
        par(new = TRUE)
        cl$object <- itr
        eval(cl)
      }
    })
  }
})


#--- '[' operator
setMethod("[", "MareyMap", function(x, i, j, ..., drop) {
  c(name = x@markerNames[i], phys = x@physicalPositions[i], gen = x@geneticDistances[i])
})


#--- '[[' operator
setMethod("[[", "MareyMap", function(x, i, j, ..., drop) {
  lst <- list(x@physicalPositions[[i]], x@geneticDistances[[i]], x@markerValidity[[i]])
  lst <- c(x@markerNames[[i]], lst)
  names(lst) <- c("name", "phys", "gen", "valid")
  lst
})


#--- '[[<-' operator
setMethod("[[<-", "MareyMap", function(x, i, j, ..., drop, value) {
  x@markerNames[[i]] <- value[["name"]]
  x@physicalPositions[[i]] <- value[["phys"]] 
  x@geneticDistances[[i]] <- value[["gen"]]
  x@markerValidity[[i]] <- value[["valid"]]
  for(itr in interpolations(x)) {
    interpolate(itr, x)
    x <- x + itr
  }
  x
})


# as MareyMap --> data.frame
setAs("MareyMap", "data.frame", function(from, to) {
  df <- data.frame(from@markerNames, from@physicalPositions, from@geneticDistances, from@markerValidity)
  df <- cbind(from@setName, from@mapName, df)
  names <- c("set", "map", "mkr", "phys", "gen", "vld")
  lapply(interpolations(from), function(itr) {        
    if(persistent(itr)) {
      df <<- cbind(df, rates(itr))
      names <<- c(names, name(itr))    
    }
  })
  names(df) <- names
  df
})


# as data.frame --> MareyMap
setAs("data.frame", "MareyMap", function(from, to) {
  object <- new("MareyMap")
  mapName(object) <- as.character(as.vector(from$map)[[1]])
  setName(object) <- as.character(as.vector(from$set)[[1]])
  markerNames(object) <- as.character(as.vector(from$mkr))
  geneticDistances(object) <- from$gen
  physicalPositions(object) <- from$phys
  if("vld" %in% names(from))
    markerValidity(object) <- as.vector(from$vld)
  else
    markerValidity(object) <- rep(TRUE, nrow(from)) 
  object
})


setMethod("removeMarker", c("MareyMap", "integer"), function(object, value) {
  if((value > length(object@markerNames)) || (value <= 0))
    return(object)
  df <- as(object, "data.frame")
  itrs <- interpolations(object)
  if (value == 1)
    object <- as(df[2:nrow(df), ], "MareyMap")
  if (value == nrow(df))
    object <- as(df[1:(nrow(df) - 1), ], "MareyMap")
  if(value > 1 && value < nrow(df))
    object <- as(rbind(df[1:(value - 1), ], df[(value + 1):nrow(df), ]), "MareyMap")
  for (itr in itrs) {
    object <- object + itr
  }
  object
})


setMethod("textFile", c("MareyMap", "character"), function(object, file) {
  df <- as(object, "data.frame")
  itrcreate <- ""
  lapply(interpolations(object), function(itr) {
    if(persistent(itr))
      itrcreate <<- paste(itrcreate, "#", setName(object), " - ", mapName(object) ," - ", name(itr), ":", createOrder(itr), "\n", sep = "")
  })
  write(itrcreate, file = as.character(file))
  write.table(t(as.vector(colnames(df))), file = as.character(file), col.names = FALSE, row.names = FALSE, append = TRUE)
  write.table(df, file = as.character(file), col.names = FALSE, row.names = FALSE, append = TRUE)
})


setMethod("query", c("MareyMap", "numeric"), function(object, pos) {
  sapply(interpolations(object), function(itr) {query(itr, pos)})
})


setAs("MareyMap", "NULL", function(from, to) {NULL})

