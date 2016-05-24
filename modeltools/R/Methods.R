
### Definition and construction method for class `ModelEnv'

setMethod("initialize", signature = "ModelEnv", 
    definition = function(.Object) {

        ### a new environment: all data are stored here 
        .Object@env <- new.env()
    
        ### extract a variable names `which' from the environment
        .Object@get <-
            function(which) get(which, envir = .Object@env, 
                                inherits = FALSE)

        ### set a variable
        .Object@set <-
            function(which, data) assign(which, data, .Object@env)

        return(.Object)
    }
)

### some utility methods for ModelEnv onjects

setMethod("show", signature = "ModelEnv", 
    definition = function(object) {

    cat("\n")
    cat("A", class(object), "with \n\n")
    n <- NULL
    if (has(object, "response")) {
        cat("  response variable(s):  ",
            colnamesnum(object@get("response")), "\n")
        n <- nrow(object@get("response"))
    }
    else if (has(object, "responseMatrix")) {
        cat("  response matrix column(s): ",
            colnamesnum(object@get("responseMatrix")), "\n")
        n <- nrow(object@get("responseMatrix"))
    }
    
    if (has(object, "input")) {
        cat("  input variable(s):     ",
            colnamesnum(object@get("input")), "\n")
        n <- nrow(object@get("input"))
    }
    else if (has(object, "designMatrix")) {
        cat("  design matrix column(s): ",
            colnamesnum(object@get("designMatrix")), "\n")
        n <- nrow(object@get("designMatrix"))
    }
    
    if (is.null(n)) 
        cat("  no observations\n")
    else
        cat("  number of observations:", n, "\n")

    if(length(object@hooks)>0){
            for(n in 1:length(object@hooks)){
                if(n==1)
                    cat("  hooks                 : ")
                else
                    cat("                          ")
                
                cat(paste(names(object@hooks)[n],"(",
                          paste(names(object@hooks[[n]]), collapse=", "),
                          ")", sep=""), "\n")
            }
        }
                    
    cat("\n")
        

})

## Utility function: return either names or number of columns
colnamesnum <- function(x)
{
    if(is.null(colnames(x)))
        return(ncol(x))
    else
        return(colnames(x))
}


setGeneric("has", function(object, which) standardGeneric("has"))

setMethod("has", signature(object = "ModelEnv", which = "character"),
    definition = function(object, which) {
        exists(which, envir = object@env, inherits = FALSE)
    }
)

setGeneric("dimension", function(object, which) standardGeneric("dimension"))

setMethod("dimension", signature(object = "ModelEnv", which = "character"),
    definition = function(object, which) {
        if (has(object, which))
            eval(parse(text = paste("dim(",which,")")) , envir = object@env)
        else
            NULL
    }
)


setGeneric("empty", function(object) standardGeneric("empty"))

setMethod("empty", signature(object = "ModelEnv"),
          definition = function(object) length(ls(object@env))==0)


###**********************************************************

setGeneric("clone", function(object, ...) standardGeneric("clone"))

## the set() method of ModelEnvFormula objects uses lexical scope on
## various bits and pieces, hence cloning currently returns only a
## ModelEnv object, which only has a trivial get method and no set
## method

setMethod("clone", signature = "ModelEnv", 
    definition = function(object, copydata = TRUE) {

    z <- new(class(object))

    if (extends(class(object), "ModelEnvFormula"))
        z@formula <- object@formula

    ### call set and get from object, however, they work in z@env
    ### <FIXME> what about parent.frame() ???
    z@set <- function(which = NULL, data = NULL, frame = parent.frame(), 
                      envir = z@env)
        object@set(which = which, data = data, frame = frame, env = envir)
    z@get <- function(which, data = NULL, frame = parent.frame(), 
                      envir = z@env)
        object@get(which = which, data = data, frame = frame, env = envir)
    ### </FIXME>

    if (copydata) {
        for (name in ls(object@env))
            assign(name, object@get(name), envir = z@env)
    }
    return(z)
})

setGeneric("subset", function(x, ...) standardGeneric("subset"))

setMethod("subset", signature = "ModelEnv",
          definition = function(x, subset, clone = TRUE, ...)
{
    MYSUBSET <- function(x, subset, ...){
        if (is(x, "matrix"))
            x[subset,,drop=FALSE]
        else
            subset(x, subset, ...)
    }

    z <- MEapply(x, MYSUBSET, clone=clone, subset=subset, ...)

    if (!clone)
        invisible(z)
    else
        return(z)
    
})

### dpp, fit and predict generics for StatModel objects

setGeneric("fit", function(model, data, ...) standardGeneric("fit"))

setMethod("fit", signature = signature(model = "StatModel", 
                                       data = "ModelEnv"),
    definition = function(model, data, ...)
        model@fit(data, ...)
)

setGeneric("dpp", function(model, ...) standardGeneric("dpp"))

setMethod("dpp", signature = "StatModel", 
    definition = function(model, ...)
        model@dpp(...)
)

### don't want to redefine stats:::predict, but ...

Predict <- function(object, ...) {
    if ("statmodel" %in% names(object)) {
        if (is(object$statmodel, "StatModel"))
            return(object$statmodel@predict(object, ...))
    }
    return(predict(object, ...))
}
