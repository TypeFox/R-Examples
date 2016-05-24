init.ShojiCatalog <- function (.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
    names(.Object@index) <- absoluteURL(names(.Object@index), .Object@self)
    return(.Object)
}
setMethod("initialize", "ShojiCatalog", init.ShojiCatalog)

is.shojiCatalog <- function (x) inherits(x, "ShojiCatalog")

getIndexSlot <- function (x, i, what=character(1)) {
    vapply(index(x), function (a) a[[i]], what, USE.NAMES=FALSE)
}

setIndexSlot <- function (x, i, value, unique=FALSE) {
    if (length(value) == 1) value <- rep(value, length(x))
    stopifnot(length(x) == length(value))

    old <- index(x)
    index(x) <- mapply(function (a, v) {
        a[[i]] <- v
        return(a)
    }, a=index(x), v=value, SIMPLIFY=FALSE)
    if (unique) {
        ## Check to see if any of the value is duplicated after updating
        newvals <- getIndexSlot(x, i) ## Assumes "character". Revisit if need unique for non-char
        dups <- duplicated(newvals)
        if (any(dups)) {
            halt("Duplicate values not permitted: ",
                serialPaste(unique(newvals[dups])))
        }
    }
    to.update <- dirtyElements(old, index(x))
    if (any(to.update)) {
        ## Make sure certain fields are [] in the JSON
        ensure <- c("subvariables")
        payload <- lapply(index(x)[to.update], function (p) {
            p <- p[i] ## Let's only PATCH the field we're changing
            these <- intersect(ensure, names(p))
            if (length(these)) p[these] <- lapply(p[these], I)
            return(p)
        })
        crPATCH(self(x), body=toJSON(payload))
    }
    return(x)
}

dirtyElements <- function (x, y) {
    !mapply(identical, x, y, USE.NAMES=FALSE, SIMPLIFY=TRUE)
}

##' @rdname catalog-extract
##' @export
setMethod("[", c("ShojiCatalog", "character"), function (x, i, ...) {
    w <- match(i, urls(x))
    if (any(is.na(w))) {
        halt("Undefined elements selected: ", serialPaste(i[is.na(w)]))
    }
    callNextMethod(x, w, value)
})
##' @rdname catalog-extract
##' @export
setMethod("[", c("ShojiCatalog", "numeric"), function (x, i, ...) {
    bad <- abs(as.integer(i)) > length(x)
    if (any(bad)) {
        halt("Subscript out of bounds: ", capture.output(dput(i[bad])))
    }
    callNextMethod(x, i, value)
})
##' @rdname catalog-extract
##' @export
setMethod("[", c("ShojiCatalog", "logical"), function (x, i, ...) {
    if (length(i) > length(x)) {
        halt("Subscript out of bounds: got ", length(i), " logicals, need ",
            length(x))
    }
    index(x) <- index(x)[i]
    return(x)
})
##' @rdname catalog-extract
##' @export
setMethod("[", c("ShojiCatalog", "ANY"), function (x, i, ...) {
    index(x) <- index(x)[i]
    return(x)
})
##' @rdname catalog-extract
##' @export
setMethod("[[", c("ShojiCatalog", "ANY"), function (x, i, ...) {
    index(x)[[i]]
})

##' @rdname catalog-extract
##' @export
setMethod("[[", c("ShojiCatalog", "character"), function (x, i, ...) {
    stopifnot(length(i) == 1L)
    w <- whichNameOrURL(x, i)
    if (is.na(w)) {
        return(NULL)
    }
    index(x)[[w]]
})

##' @rdname catalog-extract
##' @export
setMethod("$", "ShojiCatalog", function (x, name) x[[name]])

##' Length of Catalog
##' @param x a Catalog
##' @return Integer: the number of elements in the index list
##' @name catalog-length
NULL

whichNameOrURL <- function (x, i) {
    ns <- names(x)
    w <- match(i, ns)
    if (any(is.na(w))) {
        w <- match(i, urls(x))
    } else {
        ## Warn if duplicated
        dups <- i %in% ns[duplicated(ns)]
        if (any(dups)) {
            bads <- i[dups]
            msg <- ifelse(length(bads) > 1,
                " do not uniquely identify elements. Returning the first matches",
                " does not uniquely identify elements. Returning the first match")
            warning(i, msg, call.=FALSE)
        }
    }
    return(w)
}

##' @rdname catalog-length
##' @export
setMethod("length", "ShojiCatalog", function (x) length(index(x)))
setMethod("lapply", "ShojiCatalog", function (X, FUN, ...) lapply(index(X), FUN, ...))

##' Get the body of a Catalog
##'
##' The core of Catalog data is in its "index". These methods get and set that
##' slot.
##' @param x a Catalog (VariableCatalog, Subvariables, or similar object)
##' @param value For the setters, an appropriate-length list to
##' assign
##' @return Getters return the list object in the "index" slot; setters
##' return \code{x} duly modified.
##' @aliases index index<-
##' @name index
NULL

##' @rdname index
##' @export
setMethod("index", "ShojiCatalog", function (x) x@index)
##' @rdname index
##' @export
setMethod("index<-", "ShojiCatalog", function (x, value) {
    x@index <- value
    return(x)
})

##' Get the URLs contained in a Catalog or Order object
##'
##' Sometimes it is useful to extract flattened vector of URLs from more
##' complex objects for purposes like subsetting or doing set comparisons.
##'
##' @param x a Catalog, Order, or Group object
##' @return A character vector of URLs
##' @aliases urls
##' @keywords internal
##' @rdname urls
##' @export
setMethod("urls", "ShojiCatalog", function (x) names(index(x)))

##' @rdname describe-catalog
##' @export
setMethod("names", "ShojiCatalog", function (x) getIndexSlot(x, "name"))
##' @export
##' @rdname describe-catalog
setMethod("names<-", "ShojiCatalog", function (x, value) {
    setIndexSlot(x, "name", value, unique=TRUE)
})

##' @export
as.list.ShojiCatalog <- function (x, ...) lapply(names(index(x)), function (i) x[[i]])

##' Utility to get a more human-readable view of a Shoji Catalog
##'
##' @param x ShojiCatalog or subclass
##' @param keys character vector of attribute names from each catalog tuple to
##' include in the result. Default is TRUE, which means all.
##' @param rownames See \code{\link[base]{data.frame}}, the \code{row.names}
##' argument, to which this is passed in \code{data.frame}. The difference here
##' is that if \code{rownames} is explicitly set as \code{NULL}, the resulting
##' object will not have row names set. By default, row names will be the URLs
##' of the catalog tuples.
##' @param ... additional arguments passed to \code{data.frame}
##' @return a \code{data.frame} view of the catalog
##' @export
catalogToDataFrame <- function (x, keys=TRUE, rownames, ...) {
    default.rownames <- missing(rownames)
    if (default.rownames) {
        rownames <- NULL
    }
    out <- data.frame(do.call(rbind, lapply(index(x), function (a) a[keys])),
        row.names=rownames, ...)
    if (default.rownames) {
        rownames(out) <- NULL
    }
    return(out)
}
