##' @rdname toVariable
##' @export
setMethod("toVariable", "character", function (x, ...) {
    return(structure(list(values=x, type="text", ...),
        class="VariableDefinition"))
})
##' @rdname toVariable
##' @export
setMethod("toVariable", "numeric", function (x, ...) {
    return(structure(list(values=x, type="numeric", ...),
        class="VariableDefinition"))
})
##' @rdname toVariable
##' @export
setMethod("toVariable", "factor", function (x, ...) {
    nlevels <- length(levels(x))
    max.categories <- getOption("crunch.max.categories")
    if (!is.null(max.categories) && nlevels > max.categories) {
        return(toVariable(as.character(x), ...))
    }
    out <- structure(list(values=as.integer(x), type="categorical",
        categories=categoriesFromLevels(levels(x)), ...),
        class="VariableDefinition")
    return(NAToCategory(out, useNA="always"))
})
##' @rdname toVariable
##' @export
setMethod("toVariable", "Date", function (x, ...) {
    return(structure(list(values=as.character(x), type="datetime",
        resolution="D", ...),
        class="VariableDefinition"))
})
##' @rdname toVariable
##' @export
setMethod("toVariable", "POSIXt", function (x, ...) {
    return(structure(list(values=strftime(x, "%Y-%m-%dT%H:%M:%OS3"),
        type="datetime",
        resolution="ms", ...),
        class="VariableDefinition"))
})

##' @rdname toVariable
##' @export
setMethod("toVariable", "VariableDefinition", function (x, ...) {
    return(updateList(x, list(...)))
})
##' @rdname toVariable
##' @export
setMethod("toVariable", "logical", function (x, ...) {
    ## Make it categorical
    out <- structure(list(values=2L-as.integer(x), type="categorical",
        categories=categoriesFromLevels(c("True", "False")),
        ...),
        class="VariableDefinition")
    return(NAToCategory(out))
})

categoriesFromLevels <- function (x) {
    return(lapply(seq_along(x), function (i) {
        list(id=i, name=x[i], numeric_value=i, missing=FALSE)
    }))
}

NAToCategory <- function (var.metadata, useNA=c("ifany", "always")) {
    useNA <- match.arg(useNA)
    if (useNA == "always" || any(is.na(var.metadata$values))) {
        var.metadata$values[is.na(var.metadata$values)] <- -1L
        var.metadata$categories[[length(var.metadata$categories)+1]] <- .no.data
    }
    return(var.metadata)
}
