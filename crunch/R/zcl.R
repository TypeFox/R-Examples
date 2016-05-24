## Functions to make ZZ9 Command Language query objects

r2zcl <- function (x) {
    ## Convert an R vector to a value/column to be sent in a ZCL request.
    ## Called inside the zcl() method.
    v <- toVariable(x)
    attributes(v$values) <- NULL

    ## Grab any "typeof" attribute that's been attached. This is so we can
    ## assert that the type we're sending matches the type of some other
    ## variable, usually another variable in our CrunchExpr
    zztype <- attr(x, "typeof")

    ## If there is a single value, call it "value". Else it is a "column" array
    if (length(x) == 1) {
        out <- list(value=v$values)
    } else {
        out <- list(column=v$values)
    }

    ## Add type information since we have it, so that ZZ9 doesn't have to guess
    if (!is.null(zztype)) {
        out$type <- zztype
    } else {
        out$type <- list(value=list(class=v$type))
    }
    return(out)
}

## Methods to convert various objects to ZCL
setMethod("zcl", "CrunchExpr", function (x) x@expression)
setMethod("zcl", "CrunchVariable", function (x) list(variable=self(x)))
setMethod("zcl", "VariableTuple", function (x) list(variable=self(x)))
setMethod("zcl", "numeric", r2zcl)
setMethod("zcl", "character", r2zcl)
setMethod("zcl", "Date", r2zcl)
setMethod("zcl", "POSIXt", r2zcl)
setMethod("zcl", "logical", function (x) {
    x[is.na(x)] <- FALSE
    out <- list(column=I(x), type=list(class="boolean"))
    return(out)
})
setMethod("zcl", "NULL", function (x) NULL)
setOldClass("zcl")
setMethod("zcl", "zcl", function (x) x)
setMethod("zcl", "list", function (x) x) ## is this a good idea?
setMethod("zcl", "CrunchFilter", function (x) x@body$expression)

typeof <- function (x, variable) {
    ## Add ZCL metadata asserting that x is the same type as variable
    if (is.character(variable)) {
        variable <- list(class=variable)
    } else if (is.variable(variable) || inherits(variable, "VariableTuple")) {
        variable <- zfunc("typeof", variable)
    }
    attr(x, "typeof") <- variable
    return(x)
}

zfunc <- function (func, ...) {
    ## Wrapper that creates ZCL function syntax
    return(list(`function`=func, args=lapply(list(...), zcl)))
}
