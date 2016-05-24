## Generating function
FixRobModel <- function(center = ParamFamily(), neighbor = ContNeighborhood()){
    if(!is(neighbor, "UncondNeighborhood"))
        stop("'neighbor' is no unconditional neighborhood")
    if(any(neighbor@radius < 0 || neighbor@radius > 1))
        stop("neighborhood radius has to be in [0, 1]")

    FRM <- new("FixRobModel")
    FRM@center <- center
    FRM@neighbor <- neighbor
    
    return(FRM)
}

## Replace methods
setReplaceMethod("neighbor", "FixRobModel", 
    function(object, value){ 
        object@neighbor <- value 
        if(!is(value, "UncondNeighborhood"))
            stop("'value' is no unconditional neighborhood")
        if(any(value@radius < 0 || value@radius > 1))
            stop("neighborhood radius has to be in [0, 1]")
        object
    })
