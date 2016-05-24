## Defining generic "Weights" with no default
setGeneric("Weights", function(object) standardGeneric("Weights"))
## Defining generic "Data" with no default
setGeneric("Solution", function(object) standardGeneric("Solution"))
## Defining generic "DrawDowns" with no default
setGeneric("DrawDowns", function(object) standardGeneric("DrawDowns"))
## Setting generic functions for plot, update
setGeneric("plot")
setGeneric("update")
