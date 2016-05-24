##' Subset datasets and extract variables
##'
##' @param x a CrunchDataset
##' @param i As with a \code{data.frame}, there are two cases: (1) if no other
##' arguments are supplied (i.e \code{x[i]}), \code{i} provides for
##' \code{as.list} extraction: columns of the dataset rather than rows. If
##' character, identifies variables to extract based on their aliases (by
##' default: set \code{options(crunch.namekey.dataset="name")} to use variable
##' names); if numeric or logical,
##' extracts variables accordingly. Alternatively, (2) if \code{j} is specified
##' (as either \code{x[i, j]} or \code{x[i,]}), \code{i} is an object of class
##' \code{CrunchLogicalExpr} that will define a subset of rows.
##' @param j columnar extraction, as described above
##' @param name columnar extraction for \code{$}
##' @param drop logical: autmatically simplify a 1-column Dataset to a Variable?
##' Default is FALSE, and the TRUE option is in fact not implemented.
##' @param ... additional arguments
##' @return \code{[} yields a Dataset; \code{[[} and \code{$} return a Variable
##' @name dataset-extract
##' @aliases dataset-extract
NULL

##' @rdname dataset-extract
##' @export
setMethod("[", c("CrunchDataset", "ANY"), function (x, i, ..., drop=FALSE) {
    x@variables <- variables(x)[i]
    readonly(x) <- TRUE ## we don't want to overwrite the big object accidentally
    return(x)
})
##' @rdname dataset-extract
##' @export
setMethod("[", c("CrunchDataset", "character"), function (x, i, ..., drop=FALSE) {
    w <- match(i, names(x))
    if (any(is.na(w))) {
        halt("Undefined columns selected: ", serialPaste(i[is.na(w)]))
    }
    callNextMethod(x, w, ..., drop=drop)
})
##' @rdname dataset-extract
##' @export
setMethod("[", c("CrunchDataset", "missing", "ANY"), function (x, i, j, ..., drop=FALSE) {
    x[j]
})

##' @rdname dataset-extract
##' @export
setMethod("[", c("CrunchDataset", "CrunchLogicalExpr", "missing"), function (x, i, j, ..., drop=FALSE) {
    f <- activeFilter(x)
    if (length(zcl(f))) {
        i <- f & i
    }
    activeFilter(x) <- i
    return(x)
})

##' @rdname dataset-extract
##' @export
setMethod("[", c("CrunchDataset", "CrunchLogicalExpr", "ANY"), function (x, i, j, ..., drop=FALSE) {
    ## Do the filtering of rows, then cols
    x <- x[i,]
    return(x[j])
})

##' @rdname dataset-extract
##' @export
setMethod("subset", "CrunchDataset", function (x, ...) {
    x[..1,]
})

##' @rdname dataset-extract
##' @export
setMethod("[[", c("CrunchDataset", "ANY"), function (x, i, ..., drop=FALSE) {
    out <- variables(x)[[i]]
    if (!is.null(out)) {
        out <- CrunchVariable(out, filter=activeFilter(x))
    }
    return(out)
})
##' @rdname dataset-extract
##' @export
setMethod("[[", c("CrunchDataset", "character"), function (x, i, ..., drop=FALSE) {
    stopifnot(length(i) == 1)
    n <- match(i, names(x))
    if (is.na(n)) {
        ## See if the variable in question is hidden
        hvars <- hidden(x)
        hnames <- findVariables(hvars, key=namekey(x), value=TRUE)
        n <- match(i, hnames)
        if (is.na(n)) {
            return(NULL)
        } else {
            ## If so, return it with a warning
            out <- hidden(x)[[n]]
            if (!is.null(out)) {
                out <- CrunchVariable(out, filter=activeFilter(x))
            }
            warning("Variable ", i, " is hidden", call.=FALSE)
            return(out)
        }
    } else {
        return(callNextMethod(x, n, ..., drop=drop))
    }
})
##' @rdname dataset-extract
##' @export
setMethod("$", "CrunchDataset", function (x, name) x[[name]])


## Things that set

.addVariableSetter <- function (x, i, value) {
    if (i %in% names(x)) {
        ## We're not adding, we're updating.
        return(.updateValues(x, i, value))
    } else {
        if (inherits(value, "VariableDefinition")) {
            ## Just update its alias with the one we're setting
            value$alias <- i
            ## But also check to make sure it has a name, and use `i` if not
            value$name <- value$name %||% i
        } else {
            ## Create a VarDef and use `i` as name and alias
            value <- VariableDefinition(value, name=i, alias=i)
        }
        addVariables(x, value)
    }
}

.updateValues <- function (x, i, value, filter=NULL) {
    if (length(i) != 1) {
        halt("Can only update one variable at a time (for the moment)")
    }
    variable <- x[[i]]
    if (is.null(filter)) {
        variable[] <- value
    } else {
        variable[filter] <- value
    }
    return(x)
}

.updateVariableMetadata <- function (x, i, value) {
    ## Confirm that x[[i]] has the same URL as value
    v <- Filter(function (a) a[[namekey(x)]] == i,
        index(allVariables(x)))
    if (length(v) == 0) {
        ## We may have a new variable, and it's not
        ## yet in our variable catalog. Let's check.
        x <- refresh(x)
        if (!(self(value) %in% urls(allVariables(x)))) {
            halt("This variable does not belong to this dataset")
        }
        ## Update value with `i` if it is
        ## different. I.e. set the alias based on i if not otherwise
        ## specified. (setTupleSlot does the checking)
        tuple(value) <- setTupleSlot(tuple(value), namekey(x), i)
    } else if (!identical(names(v), self(value))) {
        ## x[[i]] exists but is a different variable than value
        halt("Cannot overwrite one Variable with another")
    }
    allVariables(x)[[self(value)]] <- value
    return(x)
}

##' Update a variable or variables in a dataset
##'
##' @param x a CrunchDataset
##' @param i For \code{[}, a \code{CrunchLogicalExpr}, numeric, or logical
##' vector defining a subset of the rows of \code{x}. For \code{[[}, see
##' \code{j} for the as.list column subsetting.
##' @param j if character, identifies variables to extract based on their
##' aliases (by default: set \code{options(crunch.namekey.dataset="name")} 
##' to use variable names); if numeric or
##' logical, extracts variables accordingly. Note that this is the as.list
##' extraction, columns of the dataset rather than rows.
##' @param name like \code{j} but for \code{$}
##' @param value replacement values to insert. These can be \code{crunchExpr}s
##' or R vectors of the corresponding type
##' @return \code{x}, modified.
##' @aliases dataset-update
##' @name dataset-update
NULL

##' @rdname dataset-update
##' @export
setMethod("[[<-",
    c("CrunchDataset", "character", "missing", "CrunchVariable"),
    .updateVariableMetadata)
##' @rdname dataset-update
##' @export
setMethod("[[<-",
    c("CrunchDataset", "ANY", "missing", "CrunchVariable"),
    function (x, i, value) .updateVariableMetadata(x, names(x)[i], value))
##' @rdname dataset-update
##' @export
setMethod("[[<-",
    c("CrunchDataset", "character", "missing", "ANY"),
    .addVariableSetter)
##' @rdname dataset-update
##' @export
setMethod("[[<-",
    c("CrunchDataset", "character", "missing", "CrunchLogicalExpr"),
    function (x, i, value) {
        halt("Cannot currently derive a logical variable")
    })
##' @rdname dataset-update
##' @export
setMethod("[[<-",
    c("CrunchDataset", "ANY"),
    function (x, i, value) {
        halt("Only character (name) indexing supported for [[<-")
    })
##' @rdname dataset-update
##' @export
setMethod("[[<-",
    c("CrunchDataset", "character", "missing", "NULL"),
    function (x, i, value) {
        if (!(i %in% names(x))) {
            message(dQuote(i), " is not a variable; nothing to delete by assigning NULL")
            return(x)
        }
        return(deleteVariables(x, i))
    })
##' @rdname dataset-update
##' @export
setMethod("[[<-",
    c("CrunchDataset", "ANY", "missing", "NULL"),
    function (x, i, value) deleteVariables(x, names(x)[i]))
##' @rdname dataset-update
##' @export
setMethod("$<-", c("CrunchDataset"), function (x, name, value) {
    x[[name]] <- value
    return(x)
})

##' @rdname dataset-update
##' @export
setMethod("[<-", c("CrunchDataset", "ANY", "missing", "list"),
    function (x, i, j, value) {
        ## For lapplying over variables to edit metadata
        stopifnot(length(i) == length(value),
            all(vapply(value, is.variable, logical(1))))
        for (z in seq_along(i)) {
            x[[i[z]]] <- value[[z]]
        }
        return(x)
    })

## TODO: add similar [<-.CrunchDataset, CrunchDataset/VariableCatalog

##' @rdname dataset-update
##' @export
setMethod("[<-", c("CrunchDataset", "CrunchExpr", "ANY", "ANY"),
     function (x, i, j, value) {
        if (j %in% names(x)) {
            return(.updateValues(x, j, value, filter=i))
        } else {
            halt("Cannot add variable to dataset with a row index specified")
        }
    })
