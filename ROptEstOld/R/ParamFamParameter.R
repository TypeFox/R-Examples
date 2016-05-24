## generating function
ParamFamParameter <- function(name, main = numeric(0), 
                              nuisance, trafo){
    if(missing(name))
        name <- "parameter of a parametric family of probability measures"
    if(missing(nuisance))
        nuisance <- NULL
    if(missing(trafo))
        trafo <- diag(length(main)+length(nuisance))

    dimension <- length(main) + length(nuisance)
    if(ncol(trafo) != dimension)
        stop("invalid transformation:\n", 
             "number of columns of 'trafo' not equal to ", 
             "dimension of the parameter")
    if(nrow(trafo) > dimension)
        stop("invalid transformation:\n",
             "number of rows of 'trafo' larger than ", 
             "dimension of the parameter")
    if(any(!is.finite(trafo)))
        stop("infinite or missing values in 'trafo'")
        
    PFP <- new("ParamFamParameter")
    PFP@name <- name
    PFP@main <- main
    PFP@nuisance <- nuisance
    PFP@trafo <- trafo
    
    return(PFP)
}

## access methods
setMethod("main", "ParamFamParameter", function(object) object@main)
setMethod("nuisance", "ParamFamParameter", function(object) object@nuisance)
setMethod("trafo", "ParamFamParameter", function(object) object@trafo)

## replace methods
setReplaceMethod("main", "ParamFamParameter", 
    function(object, value){ 
        object@main <- value
        dimension <- length(object@main) + length(object@nuisance)
        if(ncol(object@trafo) != dimension)
            warning("number of columns of 'trafo' not equal to dimension of the parameters\n")
        if(nrow(object@trafo) > dimension)
            warning("number of rows of 'trafo' larger than dimension of the parameters\n")
        object
    })
setReplaceMethod("nuisance", "ParamFamParameter", 
    function(object, value){ 
        object@nuisance <- value
        dimension <- length(object@main) + length(object@nuisance)
        if(ncol(object@trafo) != dimension)
            warning("number of columns of 'trafo' not equal to dimension of the parameters\n")
        if(nrow(object@trafo) > dimension)
            warning("number of rows of 'trafo' larger than dimension of the parameters\n")
        object
    })
setReplaceMethod("trafo", "ParamFamParameter", 
    function(object, value){ 
        object@trafo <- value
        if(any(!is.finite(value)))
            stop("infinite or missing values in 'value'")
        dimension <- length(object@main) + length(object@nuisance)
        if(ncol(value) != dimension)
            warning("number of columns of 'value' not equal to dimension of the parameters\n")
        if(nrow(value) > dimension)
            warning("number of rows of 'value' larger than dimension of the parameters\n")
        object
    })

## method length
setMethod("length", "ParamFamParameter", 
    function(x){ length(x@main) + length(x@nuisance) })
