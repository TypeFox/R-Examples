## generating function
RealRandVariable <- function(Map = list(function(x){1}), Domain = NULL, Range) {
    nrvalues <- length(Map)
    for(i in 1:nrvalues){
        if(!is.function(Map[[i]])) 
            stop("element ", i, " of 'Map' contains no function")
        if(length(formals(Map[[i]])) != 1)
            stop("element ", i, " of 'Map' has to be a function of one argument")
        if(names(formals(Map[[i]])) != "x")
            stop("element ", i, " of 'Map' contains a function with argument name != 'x'")
    }

    if(missing(Range)) Range <- Reals()
    if(!is(Range, "Reals"))
        stop("'Range' has to be of class 'Reals'")

    R <- new("RealRandVariable")
    R@Map <- Map
    R@Domain <- Domain
    R@Range <- Range
    
    return(R)
}

## replace method
setReplaceMethod("Range", "RealRandVariable", 
    function(object, value){ 
        object@Range <- value 
        if(!is(value, "Reals"))
            stop("'Range' of 'value' is not the Real space")
        object
    })
