## generating function
ParamFamily <- function(name, distribution = Norm(), distrSymm, 
                        main = 0, nuisance, trafo, param, 
                        props = character(0)){
    if(missing(name)) 
        name <- "parametric family of probability measures"
    if(missing(distrSymm)) distrSymm <- NoSymmetry()
    if(missing(param)) 
        param <- ParamFamParameter(name = paste("parameter of", name), 
                        main = main, nuisance = nuisance, trafo = trafo)
    PF <- new("ParamFamily")
    PF@name <- name
    PF@distribution <- distribution
    PF@distrSymm <- distrSymm
    PF@param <- param
    PF@props <- props
    
    return(PF)
}

## access methods
setMethod("param", "ParamFamily", function(object) object@param)

## wrapped access methods
setMethod("main", "ParamFamily", function(object) main(param(object)))
setMethod("nuisance", "ParamFamily", function(object) nuisance(param(object)))
setMethod("trafo", "ParamFamily", function(object) trafo(param(object)))

## replace methods
#setReplaceMethod("param", "ParamFamily", 
#    function(object, value){ object@param <- value; object })
