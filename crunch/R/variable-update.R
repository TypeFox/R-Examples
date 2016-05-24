.updateVariable <- function (variable, value, filter=NULL) {
    ## Construct a ZCL update payload, then POST it
    payload <- list(command="update",
        variables=.updatePayload(variable, value))
    payload[["filter"]] <- zcl(filter)
    dref <- datasetReference(variable)
    update_url <- paste0(dref, "table/")
    out <- crPOST(update_url, body=toJSON(payload))
    dropCache(dref)
    invisible(out)
}

.updatePayload <- function (variable, value) {
    ## Construct the "variables" key of a ZCL update payload
    if (is.Array(variable)) {
        out <- sapply(urls(subvariables(variable)), function (x) {
            zcl(typeof(value,
                structure(zfunc("typeof",
                structure(list(variable=x), class="zcl")),
                    class="zcl")))
            }, simplify=FALSE)
    } else {
        out <- structure(list(zcl(typeof(value, variable))),
            .Names=self(variable))
    }

    ## Check for missingness and replace the NAs with special values
    out <- lapply(out, function (x) {
        if ("column" %in% names(x)) {
            missings <- is.na(x$column)
            if (any(missings)) {
                x$column <- as.list(x$column)
                x$column[missings] <- rep(list(.no.data.value(type(variable))), sum(missings))
            }
        }
        return(x)
    })
    return(out)
}

.dispatchFilter <- function (f) {
    ## Given a valid R index (numeric, logical) or CrunchExp, make a ZCL (?) filter
    if (is.logical(f)) {
        ## Validate
        f <- which(f)
    }
    if (is.numeric(f)) {
        ## Validate

        f <- .seqCrunch(zfunc("row"), f - 1)
    }
    return(f)
}

##' Updating variables with expressions or values
##'
##' @param x a Variable
##' @param i a CrunchLogicalExpr or R index, optionally
##' @param j Invalid
##' @param value an R vector or a CrunchExpr with which to update
##' @return \code{x} duly modified
##' @name variable-update
##' @aliases variable-update
NULL

##' @rdname variable-update
##' @export
setMethod("[<-", c("CrunchVariable", "ANY", "missing", "ANY"), .backstopUpdate)

##' @rdname variable-update
##' @export
setMethod("[<-", c("CrunchVariable", "ANY", "missing", "NULL"),
    function (x, i, j, value) return(NULL))

.var.updater <- function (x, i, j, value) {
    if (missing(i)) i <- NULL
    out <- .updateVariable(x, value, filter=.dispatchFilter(i))
    return(x)
}

##' @rdname variable-update
##' @export
setMethod("[<-", c("TextVariable", "ANY", "missing", "character"),
    .var.updater)
##' @rdname variable-update
##' @export
setMethod("[<-", c("NumericVariable", "ANY", "missing", "numeric"),
    .var.updater)
##' @rdname variable-update
##' @export
setMethod("[<-", c("DatetimeVariable", "ANY", "missing", "Date"),
    .var.updater)
##' @rdname variable-update
##' @export
setMethod("[<-", c("DatetimeVariable", "ANY", "missing", "POSIXt"),
    .var.updater)
##' @rdname variable-update
##' @export
setMethod("[<-", c("CrunchVariable", "ANY", "missing", "CrunchExpr"),
    .var.updater)
##' @rdname variable-update
##' @export
setMethod("[<-", c("CrunchVariable", "CrunchExpr", "missing", "CrunchExpr"),
    function (x, i, j, value) {
        if (!identical(zcl(i), value@filter)) {
            halt("Cannot update a variable with a value that has a different filter")
        } else {
            callNextMethod()
        }
    })

## Set of functions to use in multiple dispatches
.categorical.update <- list(
    numeric=function (x, i, j, value) {
        if (missing(i)) i <- NULL
        if (all(c(NA, -1) %in% value)) {
            halt("Cannot have both NA and -1 when specifying category ids")
        }
        value[is.na(value)] <- -1
        invalids <- setdiff(value, ids(categories(x)))
        add.no.data <- -1 %in% invalids
        invalids <- setdiff(invalids, -1)
        if (length(invalids)) {
            plural <- length(invalids) > 1
            halt(paste0("Input value", ifelse(plural, "s ", " "),
                serialPaste(invalids), ifelse(plural, " are ", " is "),
                "not present in the category ids of variable ", dQuote(name(x))))
        }
        if (add.no.data) {
            categories(x)[[length(categories(x)) + 1]] <- Category(data=.no.data)
        }
        out <- .updateVariable(x, value, filter=.dispatchFilter(i))
        return(x)
    },
    character=function (x, i, j, value) {
        if (missing(i)) i <- NULL
        value[is.na(value)] <- "No Data"
        invalids <- setdiff(value, names(categories(x)))
        add.no.data <- "No Data" %in% invalids
        invalids <- setdiff(invalids, "No Data")
        if (length(invalids)) {
            plural <- length(invalids) > 1
            halt(paste0("Input value", ifelse(plural, "s ", " "),
                serialPaste(invalids), ifelse(plural, " are ", " is "),
                "not present in the category names of variable ",
                dQuote(name(x))))
        }
        if (add.no.data) {
            categories(x)[[length(categories(x)) + 1]] <- Category(data=.no.data)
        }
        value <- n2i(value, categories(x))
        out <- .updateVariable(x, value, filter=.dispatchFilter(i))
        return(x)
    },
    factor=function (x, i, j, value) {
        if (missing(i)) i <- NULL
        x[i] <- as.character(value) ## Inefficient, but probably fine
        return(x)
    }
)

##' @rdname variable-update
##' @export
setMethod("[<-", c("CategoricalVariable", "ANY", "missing", "numeric"),
    .categorical.update[["numeric"]])
##' @rdname variable-update
##' @export
setMethod("[<-", c("CategoricalVariable", "ANY", "missing", "character"),
    .categorical.update[["character"]])
##' @rdname variable-update
##' @export
setMethod("[<-", c("CategoricalVariable", "ANY", "missing", "factor"),
    .categorical.update[["factor"]])
##' @rdname variable-update
##' @export
setMethod("[<-", c("CategoricalArrayVariable", "ANY", "missing", "numeric"),
    .categorical.update[["numeric"]])
##' @rdname variable-update
##' @export
setMethod("[<-", c("CategoricalArrayVariable", "ANY", "missing", "character"),
    .categorical.update[["character"]])
##' @rdname variable-update
##' @export
setMethod("[<-", c("CategoricalArrayVariable", "ANY", "missing", "factor"),
    .categorical.update[["factor"]])

##' @rdname variable-update
##' @export
setMethod("[<-", c("CrunchVariable", "ANY", "missing", "logical"),
    function (x, i, j, value) {
        ## For assigning NA
        if (all(is.na(value))) {
            value <- .no.data.value(type(x), add.type=TRUE)
        } else {
            ## halt()
            .backstopUpdate(x, i, j, value)
        }

        ## Datetime not yet supported, apparently
        if (is.Datetime(x)) {
            .backstopUpdate(x, i, j, value)
        } else if (is.Categorical(x) || is.CA(x) || is.MR(x)) {
            return(.categorical.update[["numeric"]](x, i, j, value))
        }
        if (missing(i)) i <- NULL
        out <- .updateVariable(x, value, filter=.dispatchFilter(i))
        return(x)
    })

.no.data.value <- function (x, add.type=FALSE) {
    if (x %in% c("categorical", "multiple_response", "categorical_array")) {
        return(-1L)
    } else {
        out <- list(`?`=-1L)
        if (add.type) {
            out <- list(value=out, type=list(class=x))
        }
        return(out)
    }
}

##' @rdname variable-update
##' @export
setMethod("is.na<-", "CrunchVariable", function (x, value) {
    ## Temporarily kill this method until API supports correctly
    halt("is.na<- not yet supported for CrunchVariables")

    lab <- gsub('"', "", deparse(substitute(value)))
    value <- zcl(.dispatchFilter(value))
    payload <- structure(list(value), .Names=lab)
    # cat(toJSON(payload))
    out <- crPOST(x@fragments$missing_rules, body=toJSON(payload))
    return(x)
})
