setMethod("width", signature(object="Lattice"), function(object) object@width)
setMethod("Length", signature(object="Lattice"), function(object) object@Length)
setMethod("pivot", signature(object="Lattice"), function(object) object@pivot)

setReplaceMethod("width",  signature(object="Lattice"), 
          function(object, value) {object@width <- value; object})
setReplaceMethod("Length", signature(object="Lattice"), 
          function(object, value) {object@Length <- value; object})
setReplaceMethod("pivot",  signature(object="Lattice"), 
          function(object, value) {object@pivot <- value; object})

