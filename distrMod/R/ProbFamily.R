### from Matthias' thesis/ROptEst
## access methods
setMethod("name", "ProbFamily", function(object) object@name)
setMethod("distribution", "ProbFamily", function(object) object@distribution)
setMethod("distrSymm", "ProbFamily", function(object) object@distribution@Symmetry)
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

# wrapped accessors to slots of slot distribution
setMethod("r", "ProbFamily", function(object) r(distribution(object)))
setMethod("d", "ProbFamily", function(object) d(distribution(object)))
setMethod("p", "ProbFamily", function(object) p(distribution(object)))
setMethod("q", "ProbFamily", function(save = "default", status = 0, 
                              runLast = TRUE) q(distribution(save)))
           ### odd arg-list due to existing function in base package 
