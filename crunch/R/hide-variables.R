##' Hide and Unhide Variables
##' @param x a Variable or subset of a VariableCatalog to hide or unhide
##' @return (invisibly) the Variable or VariableCatalog, hidden or unhidden
##' @name hide
##' @aliases hide unhide
NULL

##' @rdname hide
##' @export
setMethod("hide", "CrunchVariable", function (x) {
    invisible(setTupleSlot(x, "discarded", TRUE))
})
##' @rdname hide
##' @export
setMethod("hide", "VariableCatalog", function (x) {
    invisible(setIndexSlot(x, "discarded", TRUE))
})

##' @rdname hide
##' @export
setMethod("unhide", "CrunchVariable", function (x) {
    invisible(setTupleSlot(x, "discarded", FALSE))
})
##' @rdname hide
##' @export
setMethod("unhide", "VariableCatalog", function (x) {
    invisible(setIndexSlot(x, "discarded", FALSE))
})

##' Hide and Unhide Variables Within a Dataset
##' @param dataset the Dataset to modify
##' @param x same as \code{dataset}, for `hiddenVariables<-`
##' @param variables names or indices of variables to (un)hide
##' @param value same as \code{variables}, for `hiddenVariables<-`
##' @param pattern optional regular expression to identify Variables to (un)hide
##' @param key the Variable attribute to \code{\link{grep}} with the
##' \code{pattern}. Default is "alias"
##' @param ... optional additional arguments to \code{grep}
##' @return (invisibly) \code{dataset} with the specified variables (un)hidden
##' @seealso \code{\link{hide}}
##' @export
hideVariables <- function (dataset, variables=NULL, pattern=NULL, key=namekey(dataset), ...) {
    var.urls <- findVariableURLs(dataset, refs=variables, pattern=pattern, key=key, ...)
    allVariables(dataset)[var.urls] <- hide(allVariables(dataset)[var.urls])
    invisible(dataset)
}

##' @rdname hideVariables
##' @export
`hiddenVariables<-` <- function (x, value) {
    if (is.character(value)) {
        value <- na.omit(match(value, names(x)))
    }
    if (length(value)) {
        return(hideVariables(x, value))
    } else {
        return(x)
    }
}

##' @rdname hideVariables
##' @export
unhideVariables <- function (dataset, variables=NULL, pattern=NULL,
                            key=namekey(dataset), ...) {
    var.urls <- findVariableURLs(hidden(dataset), refs=variables, pattern=pattern, key=key, ...)
    allVariables(dataset)[var.urls] <- unhide(allVariables(dataset)[var.urls])
    invisible(dataset)
}

##' Show the names of hidden variables within the dataset
##' @param dataset the Dataset
##' @param key the Variable attribute to return. Default is "alias"
##' @return a vector of the names of Variables marked as hidden.
##' @export
hiddenVariables <- function (dataset, key="name") {
    hv <- hidden(dataset)
    if (length(hv)) {
        return(sort(vapply(index(hv), function (x) x[[key]], character(1),
            USE.NAMES=FALSE)))
    } else {
        return(c())
    }
}

##' Delete Variables Within a Dataset
##' @param dataset the Dataset to modify
##' @param variables names or indices of variables to delete
##' @param pattern optional regular expression to identify Variables to delete
##' @param key the Variable attribute to \code{\link{grep}} with the
##' \code{pattern}. Default is "alias"
##' @param confirm logical: should the user be asked to confirm deletion.
##' Default is \code{TRUE} if in
##' an interactive session. You can avoid the confirmation prompt if you delete
##' \code{with(\link{consent})}.
##' @param ... optional additional arguments to \code{grep}
##' @return (invisibly) \code{dataset} with the specified variables deleted
##' @seealso \code{\link{hide}}
##' @export
deleteVariables <- function (dataset, variables=NULL, pattern=NULL, key=namekey(dataset), confirm=requireConsent(), ...) {
    var.urls <- findVariableURLs(dataset, refs=variables, pattern=pattern, key=key, ...)
    if (length(var.urls) == 1) {
        varnames <- names(allVariables(dataset)[var.urls])
        prompt <- paste0("Really delete ", dQuote(varnames), "?")
    } else {
        prompt <- paste0("Really delete these ", length(var.urls),
            " variables?")
    }
    if (confirm && !askForPermission(prompt)) {
        halt("Must confirm deleting variable(s)")
    }
    out <- lapply(var.urls, function (x) try(crDELETE(x)))
    ## Now, delete subvariables. When deleting an array, the DELETE request
    ## returns the subvariable URLs as a response body
    subvar.urls <- unlist(Filter(Negate(is.error), out))
    out2 <- lapply(subvar.urls, function (x) try(crDELETE(x)))
    invisible(refresh(dataset))
}

##' @rdname deleteVariables
##' @export
deleteVariable <- deleteVariables
