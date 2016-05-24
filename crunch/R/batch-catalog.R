init.BatchCatalog <- function (.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
    .Object@index <- .Object@index[order(names(.Object@index))]
    return(.Object)
}
setMethod("initialize", "BatchCatalog", init.BatchCatalog)

setMethod("imported", "BatchCatalog", function (x) {
    index(x) <- Filter(function (a) isTRUE(a$status == "imported"), index(x))
    return(x)
})

setMethod("pending", "BatchCatalog", function (x) {
    index(x) <- Filter(function (a) !isTRUE(a$status == "imported"), index(x))
    return(x)
})

##' @rdname describe-catalog
##' @export
setMethod("names", "BatchCatalog", function (x) urls(x))
