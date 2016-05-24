## generating function
RandVariable <- function(Map = list(function(x){ }), Domain = NULL, Range = NULL){
    nrvalues <- length(Map)
    for(i in 1:nrvalues){
        if(!is.function(Map[[i]])) 
            stop("element ", i, " of 'Map' contains no function")
        if(length(formals(Map[[i]])) != 1)
            stop("element ", i, " of 'Map' has to be a function of one argument")
        if(names(formals(Map[[i]])) != "x")
            stop("element ", i, " of 'Map' contains a function with argument name != 'x'")
    }

    R <- new("RandVariable")
    R@Map <- Map
    R@Domain <- Domain
    R@Range <- Range
    
    return(R)
}

## access methods
setMethod("Map", "RandVariable", function(f, ...) f@Map)
setMethod("Domain", "RandVariable", function(object) object@Domain)
setMethod("Range", "RandVariable", function(object) object@Range)

## Replace Methoden
setReplaceMethod("Map", "RandVariable", 
    function(object, value){ 
        object@Map <- value 
        nrvalues <- length(value)
        for(i in 1:nrvalues){
            if(!is.function(value[[i]])) 
                stop("'value' is no list of functions")                
            if(length(formals(value[[i]])) > 1)
                stop("'value' contains a function with more than one argument")
            if(names(formals(value[[i]])) != "x")
                stop(paste("element", i, "of 'Map' of 'value' contains a function", 
                            "with argument name != 'x'"))
        }
        object
    })
setReplaceMethod("Domain", "RandVariable", 
    function(object, value){ object@Domain <- value; object})
setReplaceMethod("Range", "RandVariable", 
    function(object, value){ object@Range <- value; object})

## method length
setMethod("length", "RandVariable", function(x){ length(x@Map) })

## compatible domains?
setMethod("compatibleDomains", signature(e1 = "RandVariable", 
                                         e2 = "RandVariable"), 
    function(e1, e2){
        if(is.null(e1@Domain)&is.null(e2@Domain))
            return(TRUE)
        if(is(e1@Domain, "EuclideanSpace")&is(e2@Domain, "EuclideanSpace"))
            if(e1@Domain@dimension == e2@Domain@dimension)
                return(TRUE)
        return(FALSE)
    })
