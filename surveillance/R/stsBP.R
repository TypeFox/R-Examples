######################################################################
# initialize-method for "stsBP" objects
######################################################################

fix.dimnamesBP <- function (x) {
    dimnames(x@ci) <- dimnames(x@lambda) <-
        c(dimnames(x@observed), list(NULL))
    x
}
    
init.stsBP <- function(.Object, ..., ci, lambda)
{
    .Object <- callNextMethod()  # use initialize,sts-method
    ## NOTE: we cannot have a validity check for the dimensions of ci and lambda
    ## in the class definition of "stsBP" since we could not easily get
    ## new("stsBP") to be a valid object. Thus, we will directly check here.

    ## check/set extra stsBP-slots
    dimObserved <- dim(.Object@observed)
    if (missing(ci)) {
        .Object@ci <- array(NA_real_, dim = c(dimObserved, 2L))
    } else {
        dimCI <- dim(.Object@ci)
        if (length(dimCI) != 3 || any(dimCI != c(dimObserved, 2L)))
            stop("dim(ci) = (", paste0(dimCI, collapse=","), ")")
    }
    if (missing(lambda)) {
        .Object@lambda <- array(NA_real_, dim = c(dimObserved, 0L))
    } else {
        dimLambda <- dim(.Object@lambda)
        if (length(dimLambda) != 3 || !identical(dimLambda[1:2], dimObserved))
            stop("dim(lambda) = (", paste0(dimLambda, collapse=","), ")")
    }

    ## fix dimnames of extra stsBP-slots
    .Object <- fix.dimnamesBP(.Object)
    
    return(.Object)
}

setMethod("initialize", "stsBP", init.stsBP)


######################################################################
# Special coerce method to account for consistent dimensions
######################################################################

setAs(from = "sts", to = "stsBP", function (from) {
    res <- new("stsBP", from,
               ci = array(NA_real_, dim = c(dim(from@observed), 2L)),
               lambda = array(NA_real_, dim = c(dim(from@observed), 0L)))
    fix.dimnamesBP(res)
})
