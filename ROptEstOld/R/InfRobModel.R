## Generating function
InfRobModel <- function(center = L2ParamFamily(), neighbor = ContNeighborhood()){
    if(!is(neighbor, "UncondNeighborhood"))
        stop("'neighbor' is no unconditional neighborhood")
    if(any(neighbor@radius < 0))
        stop("'radius' has to be in [0, Inf]")

    IRM <- new("InfRobModel")
    IRM@center <- center
    IRM@neighbor <- neighbor
    
    return(IRM)
}

## Replace methods
setReplaceMethod("neighbor", "InfRobModel", 
    function(object, value){ 
        object@neighbor <- value 
        if(!is(value, "UncondNeighborhood"))
            stop("'value' is no unconditional neighborhood")
        if(any(value@radius < 0))
            stop("'radius' has to be in [0, Inf]")
        object
    })

 
