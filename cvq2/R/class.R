#S4-class definition
setClass( "q2", representation( result="list", output="list" ) )
setMethod("show", "q2", function(object) func.output.performanceValues(object@result, object@output) )

setClass( "cvq2", contains = "q2" )

#  setMethod("[", "cvq2", function(x, i, j, ..., drop) { x@result <- x@result[i]; x })
#  setMethod("@",...) ???
