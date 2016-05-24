# Author: Babak Naimi, naimi.b@gmail.com
# Date :  March 2016
# Version 1.1
# Licence GPL v3

#-------------
setMethod('coordinates', signature(obj='sdmdata'),
          function(obj,...) {
            if (!is.null(obj@info) && !is.null(obj@info@coords)) obj@info@coords[,2:3]
          }
)

setMethod('coordinates', signature(obj='sdmModels'),
          function(obj,...) {
            if (!is.null(obj@data@info) && !is.null(obj@data@info@coords)) obj@data@info@coords[,2:3]
          }
)


if (!isGeneric("coordinates<-")) {
  setGeneric("coordinates<-", function(object,value)
    standardGeneric("coordinates<-"))
}


setReplaceMethod('coordinates', signature(object='sdmdata'), 
                 function(object,value) {
                   if (inherits(value,'matrix')) {
                     #
                   } else if (inherits(value,'character')) {
                     #
                   } else if (inherits(value,'formula')) {
                     #
                   }
                   object
                 }
)