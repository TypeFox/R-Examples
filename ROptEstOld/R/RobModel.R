## access methods
setMethod("center", "RobModel", function(object) object@center)
setMethod("neighbor", "RobModel", function(object) object@neighbor)

## replace methods
setReplaceMethod("center", "RobModel", 
    function(object, value){ object@center <- value; object })
setReplaceMethod("neighbor", "RobModel", 
    function(object, value){ object@neighbor <- value; object })
