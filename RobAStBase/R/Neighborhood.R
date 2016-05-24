## access method
setMethod("type", "Neighborhood", function(object) object@type)
setMethod("radius", "Neighborhood", function(object) object@radius)
## Replace method
setReplaceMethod("radius", "Neighborhood",
    function(object, value){object@radius <- value
                            object})
