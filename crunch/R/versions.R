init.VersionCatalog <- function (.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
    ord <- order(from8601(getIndexSlot(.Object, "last_update")),
        decreasing=TRUE)
    .Object@index <- .Object@index[ord]
    return(.Object)
}
setMethod("initialize", "VersionCatalog", init.VersionCatalog)

##' Access the saved versions of a dataset
##'
##' @param x a \code{CrunchDataset}
##' @return an object of class \code{VersionCatalog}. Supported methods on the
##' catalog include "names" and "timestamps".
##' @seealso \code{\link{saveVersion}} \code{\link{restoreVersion}}
##' @export
versions <- function (x) {
    stopifnot(is.dataset(x))
    return(VersionCatalog(crGET(shojiURL(x, "catalogs", "savepoints"))))
}

##' @rdname describe-catalog
##' @export
setMethod("names", "VersionCatalog", function (x) getIndexSlot(x, "description"))

##' @rdname describe-catalog
##' @export
setMethod("descriptions", "VersionCatalog", function (x) getIndexSlot(x, "description"))

##' @rdname describe-catalog
##' @export
setMethod("timestamps", "VersionCatalog", function (x) from8601(getIndexSlot(x, "last_update")))

##' Create a new saved version
##'
##' @param dataset a \code{CrunchDataset}
##' @param description character name to give the saved version, as in a
##' commit message. You are encouraged, though not strictly required, to give
##' versions unique descriptions.
##' @return invisibly, the URL of the newly created version
##' @seealso \code{\link{versions}} \code{\link{restoreVersion}}
##' @export
saveVersion <- function (dataset, description=NULL) {
    u <- shojiURL(dataset, "catalogs", "savepoints")
    out <- crPOST(u, body=toJSON(list(description=description)))
    invisible(dataset)
}

##' Restore a dataset to a previously saved version
##'
##' @param dataset a \code{CrunchDataset}
##' @param version either the name ("description") of the version to restore to
##' or the integer index of the version, as given by \code{versions(dataset)}
##' @return \code{dataset}, rolled back to \code{version}.
##' @seealso \code{\link{versions}} \code{\link{saveVersion}}
##' @export
restoreVersion <- function (dataset, version) {
    vcat <- versions(dataset)
    if (is.numeric(version)) {
        ## Allow passing by index
        v <- version
    } else {
        v <- which(names(vcat) == version)
        if (length(v) == 0) {
            halt(dQuote(version), " does not match any available versions")
        } else if (length(v) > 1) {
            halt(dQuote(version), " matches more than one version. Cannot restore.")
        }
    }
    revert_url <- index(vcat)[[v]]$revert
    crPOST(revert_url, body=NULL)
    invisible(refresh(dataset))
}
