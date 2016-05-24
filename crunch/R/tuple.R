##' Methods for IndexTuples
##'
##' IndexTuples are objects extracted from ShojiCatalogs. They are internally
##' used.
##'
##' @param x a Tuple
##' @param name a Tuple slot to get or set
##' @param i In [[, a Tuple slot to get
##' @param ... additional arguments to [[, ignored
##' @param value What to set in a given slot
##' @param confirm For \code{delete}, whether confirmation is required. See
##' \code{\link{delete}}.
##' @name tuple-methods
##' @aliases entity
NULL

##' @rdname tuple-methods
##' @export
setMethod("refresh", "IndexTuple", function (x) {
    dropCache(x@index_url)
    catalog <- ShojiCatalog(crGET(x@index_url))
    tup <- catalog[[x@entity_url]]
    if (is.null(tup)) {
        ## Get the object type from the (sub)class name
        cls <- sub("Tuple$", "", class(x))
        if (cls == "Index") cls <- "Object"
        halt(cls, " not found. It may have been deleted.")
    }
    x@body <- tup
    return(x)
})

##' @rdname tuple-methods
##' @export
setMethod("$", "IndexTuple", function (x, name) x@body[[name]])
##' @rdname tuple-methods
##' @export
setMethod("$<-", "IndexTuple", function (x, name, value) {
    x@body[[name]] <- value
    return(x)
})
##' @rdname tuple-methods
##' @export
setMethod("[[", "IndexTuple", function (x, i) x@body[[i]])
##' @rdname tuple-methods
##' @export
setMethod("[[<-", "IndexTuple", function (x, i, value) {
    x@body[[i]] <- value
    return(x)
})

setTupleSlot <- function (x, name, value) {
    if (!inherits(x, "IndexTuple")) {
        tuple(x) <- setTupleSlot(tuple(x), name, value)
    } else if (!identical(x[[name]], value)) {
        ## Skip updating if not modified
        x[[name]] <- value
        ## NB: no readonly mode. implement later if needed.
        payload <- toJSON(structure(list(structure(list(value), .Names=name)),
            .Names=x@entity_url))
        crPATCH(x@index_url, body=payload)
    }
    invisible(x)
}

##' @rdname tuple-methods
##' @export
setMethod("self", "IndexTuple", function (x) x@entity_url)

##' @rdname tuple-methods
##' @export
setMethod("entity", "VariableTuple", function (x) {
    return(VariableEntity(crGET(x@entity_url)))
})

##' @rdname tuple-methods
##' @export
setMethod("entity", "CrunchVariable", function (x) {
    return(VariableEntity(crGET(self(x))))
})
##' @rdname tuple-methods
##' @export
setMethod("entity", "DatasetTuple", function (x) {
    return(as.dataset(crGET(x@entity_url), tuple=x))
})

##' @rdname tuple-methods
##' @export
setMethod("delete", "IndexTuple", function (x, ...) {
    crDELETE(x@entity_url, drop=dropCache(x@index_url))
})
##' @rdname tuple-methods
##' @export
setMethod("delete", "DatasetTuple", function (x, confirm=requireConsent(), ...) {
    prompt <- paste0("Really delete dataset ", dQuote(name(x)), "?")
    if (confirm && !askForPermission(prompt)) {
        halt("Must confirm deleting dataset")
    }
    out <- callNextMethod()
    updateDatasetList()
    invisible(out)
})

##' @rdname tuple-methods
##' @export
setMethod("name", "IndexTuple", function (x) x$name)

##' @rdname tuple-methods
##' @export
setMethod("type", "IndexTuple", function (x) x$type)
