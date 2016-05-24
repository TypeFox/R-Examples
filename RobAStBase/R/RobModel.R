## access methods
setMethod("center", "RobModel", function(object) object@center)
setMethod("neighbor", "RobModel", function(object) object@neighbor)
setMethod("trafo", signature(object="RobModel",param="missing"),
           function(object) trafo(object@center))

## replace methods
setReplaceMethod("center", "RobModel",
    function(object, value){ object@center <- value; object })
setReplaceMethod("neighbor", "RobModel",
    function(object, value){ object@neighbor <- value; object })
setReplaceMethod("trafo", "RobModel",
    function(object, value){ center <- object@center
                             trafo(center) <- value
                             object@center <- center
                             object})
