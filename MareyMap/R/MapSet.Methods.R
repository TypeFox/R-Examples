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


	
#--- constructor
setMethod("MapSet", "character", function(set_name) {
	object <- new("MapSet")
	object@setName <- set_name
	object
})


#--- accessors
setMethod("setName", "MapSet", function(object) object@setName)


#--- replace methods
setReplaceMethod("setName", c("MapSet", "character"), function(object, value) {
	object@setName <- value
	lapply(object@maps, function(map) {
		setName(map) <- value
	})
	object
})


#--- [[
setMethod("[[", "MapSet", function(x, i, j = NULL, ...) {
	x@maps[[i, ...]]
})


#--- [[<-
setReplaceMethod("[[", c(x = "MapSet", value = "MareyMap"), function(x, i, j = NULL, ..., value) {
	x@maps[[i, ...]] <- value
	x
})


#--- [
setMethod("[", "MapSet", function(x, i, j = NULL, ..., drop = T) {
	x@maps[i, j, ..., drop]
})


#--- $
setMethod("$", "MapSet", function(x, name) x@maps$name)


#--- +
setMethod("+", c("MapSet", "MareyMap"), function(e1, e2) {
	if(setName(e2) != setName(e1)) {
		#--- TODO: add an error message here
		print("error in '+' method for MapSet")
	}
	#if(!is.null(e1@maps[[mapName(e2)]])){
	#	#--- TODO: add an error message here
	#	return(e1)
	#}
	e1@maps[[mapName(e2)]] <- e2
	e1
})


#--- - 
setMethod("-", c("MapSet", "character"), function(e1, e2) {
  maps <- e1@maps 
  maps[[e2]] <- NULL
  e1@maps <- maps
  e1
})


#--- length
setMethod("length", "MapSet", function(x) {
  length(x@maps)
})


#--- mapNames
setMethod("mapNames", "MapSet", function(object) {
	names(object@maps)
})


setAs("MapSet", "data.frame", function(from, to) {
	itrnames <- vector()
	lapply(from@maps, function(map) {
		lapply(interpolations(map), function(itr) {
			if(persistent(itr)) {
				itrnames <<- c(itrnames, name(itr))
			}
		})
	})
	data <- NULL
	if(length(itrnames) > 0)
		itrnames <- unique(itrnames)
	lapply(from@maps, function(map) {
		df <- as(map, "data.frame")
		if(length(itrnames) > 0) {
			for(iname in itrnames) {
				if(! iname %in% names(interpolations(map))) {
					cnam <- names(df)
					df <- cbind(df, rep(NA, length(physicalPositions(map))))
					names(df) <- c(cnam, iname)
				}
			}
		}
		if(is.null(data)) {
			data <<- df
		} else {
			data <<- rbind(data, df)
		}
	})
	data
})


setMethod("textFile", c("MapSet", "character"), function(object, file) {
	df <- as(object, "data.frame")
 	itrcreate <- ""
	lapply(object@maps, function(map) {
    lapply(interpolations(map), function(itr) { 
  	  if(persistent(itr))
	    	itrcreate <<- paste(itrcreate, "#", setName(map), " - ", mapName(map) ," - ", name(itr), ":", createOrder(itr), "\n", sep = "")
    })
  })
	write(itrcreate, file = as.character(file))
  write.table(t(as.vector(colnames(df))), file = as.character(file), col.names = FALSE, row.names = FALSE, append = TRUE)
	write.table(df, file = as.character(file), col.names = FALSE, row.names = FALSE, append = TRUE)
})
	
