init.VariableCatalog <- function (.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
    h_url <- .Object@orders$hier
    if (!is.null(h_url)) {
        o <- crGET(h_url, query=list(relative="on"))
        .Object@order <- VariableOrder(o)
    }
    return(.Object)
}
setMethod("initialize", "VariableCatalog", init.VariableCatalog)

.discardedTuple <- function (x) isTRUE(x[["discarded"]])

setMethod("active", "VariableCatalog", function (x) {
    index(x) <- Filter(Negate(.discardedTuple),
        index(x)[intersect(urls(ordering(x)), urls(x))])
    return(x)
})

setMethod("hidden", "VariableCatalog", function (x) {
    index(x) <- Filter(.discardedTuple, index(x))
    return(x)
})

##' @rdname catalog-extract
##' @export
setMethod("[[", c("VariableCatalog", "character"), function (x, i, ...) {
    tup <- index(x)[[i]]
    if (!is.list(tup)) {
        halt("Subscript out of bounds: ", i)
    }
    return(VariableTuple(index_url=self(x), entity_url=i, body=tup))
})
##' @rdname catalog-extract
##' @export
setMethod("[[", c("VariableCatalog", "ANY"), function (x, i, ...) {
    tup <- index(x)[[i]]
    if (!is.list(tup)) {
        halt("Subscript out of bounds: ", i)
    }
    return(VariableTuple(index_url=self(x), entity_url=urls(x)[i],
        body=index(x)[[i]]))
})
##' @rdname catalog-extract
##' @export
setMethod("[[<-", c("VariableCatalog", "character", "missing", "VariableTuple"),
    function (x, i, j, value) {
        index(x)[[i]] <- value@body
        return(x)
    })
##' @rdname catalog-extract
##' @export
setMethod("[[<-", c("VariableCatalog", "character", "missing", "CrunchVariable"),
    function (x, i, j, value) {
        stopifnot(i == self(value))
        x[[i]] <- tuple(value)
        return(x)
    })
##' @rdname catalog-extract
##' @export
setMethod("[", c("VariableCatalog", "VariableOrder"), function (x, i, ...) {
    index(x) <- index(x)[urls(i)]
    return(x)
})
##' @rdname catalog-extract
##' @export
setMethod("[", c("VariableCatalog", "VariableGroup"), function (x, i, ...) {
    index(x) <- index(x)[urls(i)]
    return(x)
})
##' @rdname catalog-extract
##' @export
setMethod("[<-", c("VariableCatalog", "character", "missing", "VariableCatalog"),
    function (x, i, j, value) {
        ## Validate!
        index(x)[i] <- index(value)[i]
        ## No save, I don't think. PATCH outside this fn?
        return(x)
    })
##' @rdname catalog-extract
##' @export
setMethod("[<-", c("VariableCatalog", "VariableOrder", "missing", "VariableCatalog"),
    function (x, i, j, value) {
        i <- urls(i)
        callNextMethod(x, i, value=value)
    })
##' @rdname catalog-extract
##' @export
setMethod("[<-", c("VariableCatalog", "VariableGroup", "missing", "VariableCatalog"),
    function (x, i, j, value) {
        i <- urls(i)
        callNextMethod(x, i, value=value)
    })

##' @export
##' @rdname describe-catalog
setMethod("aliases", "VariableCatalog", function (x) getIndexSlot(x, "alias"))
##' @export
##' @rdname describe-catalog
setMethod("aliases<-", "VariableCatalog", function (x, value) {
    setIndexSlot(x, "alias", value)
})

##' @export
##' @rdname describe-catalog
setMethod("descriptions", "VariableCatalog",
    function (x) getIndexSlot(x, "description"))
##' @export
##' @rdname describe-catalog
setMethod("descriptions<-", "VariableCatalog", function (x, value) {
    setIndexSlot(x, "description", value)
})

##' @export
##' @rdname describe-catalog
setMethod("types", "VariableCatalog", function (x) getIndexSlot(x, "type"))

## No setter for types<-
