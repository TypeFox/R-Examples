## access methods
setMethod("name", "ProbFamily", function(object) object@name)
setMethod("distribution", "ProbFamily", function(object) object@distribution)
setMethod("distrSymm", "ProbFamily", function(object) object@distrSymm)
setMethod("props", "ProbFamily", function(object) object@props)

## replace methods
setReplaceMethod("name", "ProbFamily", 
    function(object, value){ object@name <- value; object })
#setReplaceMethod("distribution", "ProbFamily", 
#    function(object, value){ object@distribution <- value; object })
setReplaceMethod("props", "ProbFamily", 
    function(object, value){ object@props <- value; object })

## add property
setMethod("addProp<-", "ProbFamily", 
    function(object, value){ props(object) <- c(props(object), value); object })
