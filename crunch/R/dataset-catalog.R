init.DatasetCatalog <- function (.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
    .Object@index <- .Object@index[order(selectFrom("name", .Object@index))]
    return(.Object)
}
setMethod("initialize", "DatasetCatalog", init.DatasetCatalog)

setMethod("active", "DatasetCatalog", function (x) {
    index(x) <- Filter(function (a) !isTRUE(a$archived), index(x))
    ## Not implementing user check now.
    # me <- userURL()
    # if (!is.null(me)) {
    #     index(x) <- Filter(function (a) a$owner_id == me, index(x))
    # }
    return(x)
})

setMethod("archived", "DatasetCatalog", function (x) {
    index(x) <- Filter(function (a) isTRUE(a$archived), index(x))
    ## Not implementing user check now.
    # me <- userURL()
    # if (!is.null(me)) {
    #     index(x) <- Filter(function (a) a$owner_id == me, index(x))
    # }
    return(x)
})

##' Extract and modify subsets of Catalog-type objects
##'
##' @param x a Catalog object
##' @param i which catalog elements to extract
##' @param name for \code{$}, the same as \code{i} for \code{[[}
##' @param j Invalid
##' @param drop Invalid
##' @param ... additional arguments
##' @param value For updating, an object of the appropriate class and size to
##' insert
##' @return A subset of \code{x} if extracting, otherwise \code{x} duly modified
##' @name catalog-extract
NULL

##' @rdname catalog-extract
##' @export
setMethod("[[", c("DatasetCatalog", "character"), function (x, i, ...) {
    w <- whichNameOrURL(x, i)
    x[[w]]
})
##' @rdname catalog-extract
##' @export
setMethod("[[", c("DatasetCatalog", "ANY"), function (x, i, ...) {
    b <- callNextMethod(x, i, ...)
    if (is.null(b)) return(NULL)
    DatasetTuple(index_url=self(x), entity_url=urls(x)[i],
        body=b)
})

##' Get and set names, aliases on Catalog-type objects
##'
##' These methods let you get and set names and aliases for variables in a
##' Dataset's catalog, or within \code{\link{Subvariables}} in an array
##' variable. They work like the base R names methods.
##'
##' Note that the \code{names} method on a Dataset returns the aliases of its
##' variables by default. This is controlled by
##' \code{getOption("crunch.namekey.dataset")}, which is "alias" by default.
##' Set \code{options(crunch.namekey.dataset="name")} if you wish to use
##' variable names. See the vignette on variables for more information.
##'
##' @param x a VariableCatalog, Subvariables, or similar object
##' @param value For the setters, an appropriate-length character vector to
##' assign
##' @return Getters return the character object in the specified slot; setters
##' return \code{x} duly modified.
##' @aliases describe-catalog aliases aliases<- descriptions descriptions<- types emails timestamps
##' @seealso \code{\link{Subvariables}} \code{\link{Categories}} \code{\link[base]{names}} \code{vignette("variables", package="crunch")}
##' @name describe-catalog
NULL
