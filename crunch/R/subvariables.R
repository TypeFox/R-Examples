##' Subvariables in Array Variables
##'
##' Multiple-response and categorical-array variables contain a set of
##' subvariables within them. The Subvariables class encapsulates them.
##'
##' Subvariables can be accessed from array variables (including multiple
##' response) with the \code{subvariables} method. They can be assigned back
##' with the \code{subvariables<-} setter, but there are limitations to what
##' is supported. Specifically, you can reorder subvariables, but you cannot
##' add or remove subvariables by \code{subvariables<-} assignment. See
##' \code{\link{deleteSubvariable}} to remove subvariables from an array.
##'
##' Subvariables have a \code{names} attribute that can be accessed, showing
##' the display names of the subvariables. These can be set with the
##' \code{names<-} method.
##'
##' Finally, subvariables can be accessed as regular (categorical) variables
##' with the \code{$} and \code{[[} extract methods.
##'
##' See the vignette on array variables for further details and examples.
##'
##' @param x A Variable or Subvariables object
##' @param value For the setters, the appropriate values to set
##'
##' @name Subvariables
##' @aliases Subvariables subvariables subvariables<-
##' @seealso \code{\link{subvars-extract}} \code{\link{describe-catalog}} \code{\link{deleteSubvariable}} \code{vignette("array-variables", package="crunch")}
NULL

##' @rdname Subvariables
##' @export
setMethod("subvariables", "CategoricalArrayVariable", function (x) {
    tup <- tuple(x)
    catalog_url <- absoluteURL(tup$subvariables_catalog, base=tup@index_url)
    vars <- VariableCatalog(crGET(catalog_url))
    out <- Subvariables(vars[subvariables(tup)])
    activeFilter(out) <- activeFilter(x)
    return(out)
})

##' @rdname Subvariables
##' @export
setMethod("subvariables", "VariableTuple", function (x) {
    ## Return subvariable *urls* from a Tuple, properly formatted and absolute
    return(absoluteURL(unlist(x$subvariables), base=x@index_url))
})

##' @rdname Subvariables
##' @export
setMethod("subvariables<-", c("CategoricalArrayVariable", "ANY"),
    function (x, value) {
        halt("Can only assign an object of class Subvariables")
    })
##' @rdname Subvariables
##' @export
setMethod("subvariables<-", c("CategoricalArrayVariable", "Subvariables"),
    function (x, value) {
        old <- subvariables(tuple(x))
        new <- urls(value)
        if (!setequal(old, new)) {
            halt("Can only reorder, not change, subvariables")
        }
        new <- I(new)
        if (!is.readonly(x) && !identical(new, I(old))) {
            body <- list(subvariables=new)
            payload <- toJSON(body)
            crPATCH(self(x), body=payload)
            dropCache(absoluteURL("../../cube/", self(x)))
            tuple(x)$subvariables <- new
        }
        return(x)
    })

##' @rdname describe-catalog
##' @export
setMethod("aliases", "Subvariables", function (x) {
    vapply(index(x), function (a) a$alias, character(1), USE.NAMES=FALSE)
})

##' @rdname describe-catalog
##' @export
setMethod("aliases<-", "Subvariables", function (x, value) {
    stopifnot(is.character(value), length(x) == length(value),
        !any(duplicated(value)))
    index(x) <- mapply(function (tuple, val) {
            tuple[["alias"]] <- val
            return(tuple)
        }, tuple=index(x), val=value, SIMPLIFY=FALSE, USE.NAMES=TRUE)
    crPATCH(self(x), body=toJSON(index(x)))
    return(x)
})

##' Extract and modify subsets of subvariables
##'
##' @param x Subvariables or an array Variable (which contains subvariables)
##' @param i which subvariables to extract
##' @param name For \code{$}, the name (not alias) of the subvariable to
##' extract
##' @param j Invalid
##' @param drop Invalid
##' @param ... additional arguments
##' @param value For updating, a CrunchExpr
##' @return A subset of \code{x} if extracting, otherwise \code{x} duly modified
##' @name subvars-extract
NULL

##' @rdname subvars-extract
##' @export
setMethod("[[", c("Subvariables", "character"), function (x, i, ...) {
    i <- match(i, names(x))
    if (is.na(i)) return(NULL)
    return(x[[i, ...]])
})
##' @rdname subvars-extract
##' @export
setMethod("[[", c("Subvariables", "ANY"), function (x, i, ...) {
    out <- callNextMethod(x, i, ...)
    if (!is.null(out)) {
        out <- CrunchVariable(out, filter=activeFilter(x))
    }
    return(out)
})

##' @rdname subvars-extract
##' @export
setMethod("[", c("Subvariables", "character"), function (x, i, ...) {
    w <- match(i, names(x))
    if (any(is.na(w))) {
        halt("Undefined subvariables selected: ", serialPaste(i[is.na(w)]))
    }
    return(x[w])
})

##' Delete subvariables from an array
##'
##' This function conceals the dirty work in making this happen. The array
##' gets unbound, the subvariables deleted, and then the remaining subvariable
##' are rebound into a new array.
##' @param variable the array variable
##' @param to.delete subvariable names to delete
##' @return a new version of variable without the indicated subvariables
##' @export
deleteSubvariables <- function (variable, to.delete) {
    ## Store some metadata up front
    payload <- copyVariableReferences(variable)
    subvars <- subvariables(variable)
    subvar.urls <- urls(subvars)
    subvar.names <- names(subvars)

    ## Identify subvariable URLs
    delete.these <- findVariableURLs(subvariables(variable), to.delete,
        key="name")
    ## Unbind
    all.subvar.urls <- unlist(unbind(variable))

    ## Delete
    dels <- lapply(delete.these, function (x) try(crDELETE(x)))

    ## Setdiff those deleted from those returned from unbind
    payload$subvariables <- I(setdiff(all.subvar.urls, delete.these))
    class(payload) <- "VariableDefinition"

    ## Rebind
    new_url <- POSTNewVariable(variableCatalogURL(variable), payload)

    ## Prune subvariable name prefix, or otherwise reset the names
    subvars <- Subvariables(crGET(absoluteURL("subvariables/", new_url)))
    names(subvars) <- subvar.names[match(urls(subvars), subvar.urls)]

    ## What to return? This function is kind of a hack.
    invisible(new_url)
}

##' @rdname deleteSubvariables
##' @export
deleteSubvariable <- deleteSubvariables

##' @rdname subvars-extract
##' @export
setMethod("[[<-",
    c("Subvariables", "character", "missing", "CrunchVariable"),
    function (x, i, value) {
        i <- match(i, names(x))
        if (is.na(i)) {
            ## Maybe we changed the name and that's what we're assigning back.
            ## Check URLs instead.
            i <- match(self(value), urls(x))
        }
        if (is.na(i)) {
            halt("subscript out of bounds")
        }
        x[[i]] <- value ## "callNextMethod"
        return(x)
    })
##' @rdname subvars-extract
##' @export
setMethod("[[<-",
    c("Subvariables", "ANY", "missing", "CrunchVariable"),
    function (x, i, value) {
        if (self(value) != urls(x)[i]) {
            halt("Cannot add or remove subvariables")
        }
        index(x)[[self(value)]] <- tuple(value)@body
        return(x)
    })
##' @rdname subvars-extract
##' @export
setMethod("[[<-",
    c("Subvariables", "ANY", "missing", "NULL"),
    function (x, i, value) {
        halt("Cannot add or remove subvariables")
    })
##' @rdname subvars-extract
##' @export
setMethod("[[<-",
    c("Subvariables", "ANY", "missing", "ANY"),
    function (x, i, value) {
        halt("Can only assign Variables into an object of class Subvariables")
    })
##' @rdname subvars-extract
##' @export
setMethod("$<-", c("Subvariables"), function (x, name, value) {
    x[[name]] <- value
    return(x)
})
##' @rdname subvars-extract
##' @export
setMethod("[<-", c("Subvariables", "character", "missing", "Subvariables"),
    function (x, i, value) {
        w <- match(i, names(x))
        if (any(is.na(w))) {
            halt("Undefined subvariables selected: ", serialPaste(i[is.na(w)]))
        }
        callNextMethod(x, w, value)
    })
##' @rdname subvars-extract
##' @export
setMethod("[<-", c("Subvariables", "ANY", "missing", "Subvariables"),
    function (x, i, value) {
        inbound <- vapply(value, function (a) self(a), character(1))
        if (!all(inbound %in% urls(x)[i])) {
            halt("Cannot add or remove subvariables")
        }
        index(x)[i] <- index(value)
        names(index(x))[i] <- inbound
        return(x)
    })
##' @rdname subvars-extract
##' @export
setMethod("[<-", c("Subvariables", "ANY", "missing", "ANY"),
    function (x, i, value) {
        halt("Can only assign Variables into an object of class Subvariables")
    })

##' @export
as.list.Subvariables <- function (x, ...) lapply(names(x), function (i) x[[i]])


##' @rdname describe-catalog
##' @export
setMethod("names", "CategoricalArrayVariable", function (x) {
    findVariables(subvariables(x), key=namekey(x), value=TRUE)
})

##' @rdname subvars-extract
##' @export
setMethod("[", c("CategoricalArrayVariable", "character"), function (x, i, ...) {
    w <- match(i, names(x))
    if (any(is.na(w))) {
        halt("Undefined subvariables selected: ", serialPaste(i[is.na(w)]))
    }
    return(subvariables(x)[w, ...])
})
##' @rdname subvars-extract
##' @export
setMethod("[[", c("CategoricalArrayVariable", "ANY"), function (x, i, ...) {
    return(subvariables(x)[[i, ...]])
})
##' @rdname subvars-extract
##' @export
setMethod("[[", c("CategoricalArrayVariable", "character"), function (x, i, ...) {
    i <- match(i, names(x))
    if (is.na(i)) return(NULL)
    return(x[[i, ...]])
})

##' @rdname subvars-extract
##' @export
setMethod("$", "CategoricalArrayVariable", function (x, name) x[[name]])


##' @rdname subvars-extract
##' @export
setMethod("[[<-",
    c("CategoricalArrayVariable", "ANY", "missing", "ANY"),
    function (x, i, value) {
        subvariables(x)[[i]] <- value
        return(x)
    })
##' @rdname subvars-extract
##' @export
setMethod("[[<-",
    c("CategoricalArrayVariable", "character", "missing", "ANY"),
    function (x, i, value) {
        i <- match(i, names(x))
        if (is.na(i)) {
            ## Maybe we changed the name and that's what we're assigning back.
            ## Check URLs instead.
            i <- match(self(value), urls(subvariables(x)))
        }
        if (is.na(i)) {
            halt("subscript out of bounds")
        }
        subvariables(x)[[i]] <- value ## "callNextMethod"
        return(x)
    })
##' @rdname subvars-extract
##' @export
setMethod("$<-", c("CategoricalArrayVariable"), function (x, name, value) {
    x[[name]] <- value
    return(x)
})

findParent <- function (subvar, dataset) {
    ## Utility to find the array parent, given a subvariable and its dataset
    if (is.variable(subvar)) subvar <- self(subvar)

    allvars <- index(allVariables(dataset))
    parent <- Filter(function (x) subvar %in% x$subvariables %||% c(), allvars)
    return(names(parent))
}
