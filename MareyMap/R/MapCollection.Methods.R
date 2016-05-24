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
setMethod("MapCollection", "missing", function() {
  object <- new("MapCollection")
  object
})


#--- [[
setMethod("[[", "MapCollection", function(x, i, j, ...) {
  if(missing(j))
    x@sets[[i, ...]]
  else
    x@sets[[i, ...]][[j, ...]]
})


#--- [[<-
setReplaceMethod("[[", c(x = "MapCollection", value = "MapSet"), function(x, i, j, ..., value) {
  x@sets[[i, ...]] <- value
  x
})


setReplaceMethod("[[", c(x = "MapCollection", value = "MareyMap"), function(x, i, j, ..., value) {
  set <- x@sets[[i, ...]]
  set[[j, ...]] <- value
  x[[i, ...]] <- set
  x
})


#--- [
setMethod("[", "MapCollection", function(x, i, j, ..., drop = T) {
  if(missing(j))
    x@sets[i, ..., drop]
  else
    x@sets[i, ..., drop][j, ..., drop]
})

  
#--- $
setMethod("$", "MapCollection", function(x, name) {
  x@sets$name
})


#--- +
setMethod("+", c("MapCollection", "MapCollection"), function(e1, e2) {
  lapply(e2@sets, function(sets) {
    e1 <<- e1 + sets
  })
  e1
})


#--- +
setMethod("+", c("MapCollection", "MareyMap"), function(e1, e2) {
  if(is.null(e1@sets[[setName(e2)]]))
    e1@sets[[setName(e2)]] <- MapSet(setName(e2))
  e1@sets[[setName(e2)]] <- e1@sets[[setName(e2)]] + e2
  e1
})


#--- +
setMethod("+", c("MapCollection", "MapSet"), function(e1, e2) {
  e1@sets[[setName(e2)]] <- e2
  e1
})


#--- -
setMethod("-", c("MapCollection", "character"), function(e1, e2) {
  sets <- e1@sets
  sets[[e2]] <- NULL
  e1@sets <- sets
  e1
})


#--- length
setMethod("length", "MapCollection", function(x) {
  sum <- 0
  lapply(x@sets, function(set) sum <<- sum + length(set))
  sum
})


#--- setName
setMethod("setNames", "MapCollection", function(object) {
    names(object@sets)
})


setAs("MapCollection", "data.frame", function(from, to) {
  itrnames <- vector()
  # getting the list of every interpolation names
  lapply(from@sets, function(set) {
    lapply(set@maps, function(map) {
      lapply(interpolations(map),function(itr) {
        if(persistent(itr)) {
          itrnames <<- c(itrnames, name(itr))
        }
      })
    })
  })
  # removing duplicates
  if(length(itrnames) > 0)
    itrnames <- unique(itrnames)
  data <- NULL
  lapply(from@sets,function(set) {
    df <- as(set, "data.frame")
    itrs <- NULL
    if(length(names(df)) > 6)
      itrs <- names(df)[7:length(names(df))]
    if(length(itrnames) > 0) {
      for(iname in itrnames) {
        if(! iname %in% itrs) {
          cnam <- names(df)
          df <- cbind(df, rep(NA, length(df$phys)))
          names(df) <- c(cnam, iname)
        }
      }
    }
    if(is.null(data))
      data <<- df
    else
      data <<- rbind(data, df)
  })
  data
})


setMethod("textFile", c("MapCollection", "character"), function(object, file) {
  df <- as(object, "data.frame")
  itrcreate <- ""
  lapply(object@sets, function(set) {
  	lapply(set@maps, function(map) {
      lapply(interpolations(map), function(itr) { 
        if(persistent(itr))
  		    itrcreate <<- paste(itrcreate, "#", setName(map), " - ", mapName(map) ," - ", name(itr), ":", createOrder(itr), "\n", sep = "")
      })
    })
	})
  write(itrcreate, file = as.character(file))
  write.table(t(as.vector(colnames(df))), file = as.character(file), col.names = FALSE, row.names = FALSE, append = TRUE)
  write.table(df, file = as.character(file), col.names = FALSE, row.names = FALSE, append = TRUE)
})

