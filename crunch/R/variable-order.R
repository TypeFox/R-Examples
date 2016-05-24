init.VariableOrder <- function (.Object, ..., duplicates=FALSE) {
    .Object <- callNextMethod(.Object, ...)
    dots <- list(...)
    if (length(dots) && !is.shoji(dots[[1]])) {
        .Object@graph <- .initEntities(dots, url.base=NULL)
    } else {
        .Object@graph <- .initEntities(.Object@graph, url.base=.Object@self)
    }
    duplicates(.Object) <- duplicates
    return(.Object)
}
setMethod("initialize", "VariableOrder", init.VariableOrder)

.initEntities <- function (x, url.base=NULL) {
    ## Sanitize the inputs in VariableGroup construction/updating
    ## Result should be a list, each element being either a URL (character)
    ## or VariableGroup

    ## Valid inputs:
    ## 1) A dataset: take all urls
    ## 2) A character vector of urls
    ## 3) A list of
    ## a) variables: take self
    ## b) mixed character and VariableGroups
    ## c) mixed character and lists that should be VariableGroups (fromJSON)
    if (is.dataset(x)) {
        return(.initEntities(urls(allVariables(x)), url.base=url.base))
    }
    if (is.character(x)) {
        return(.initEntities(as.list(x), url.base=url.base))
    }
    if (is.list(x)) {
        ## Init raw (fromJSON) groups
        raw.groups <- vapply(x, is.list, logical(1))
        x[raw.groups] <- lapply(x[raw.groups],
            function (a) VariableGroup(group=names(a), entities=a[[1]],
                url.base=url.base)) ## The new shoji:order structure
        ## Get self if any are Variables
        vars <- vapply(x, is.variable, logical(1))
        x[vars] <- lapply(x[vars], self)
        ## Now everything should be valid
        nested.groups <- vapply(x,
            function (a) inherits(a, "VariableGroup"),
            logical(1))
        string.urls <- vapply(x,
            function (a) is.character(a) && length(a) == 1,
            logical(1))
        stopifnot(all(string.urls | nested.groups))

        ## Absolutize if needed
        if (!is.null(url.base)) {
            x[string.urls] <- lapply(x[string.urls], absoluteURL,
                base=url.base)
        }
        ## Make sure there are no names on the list--will throw off toJSON
        names(x) <- NULL
        return(x)
    }
    halt(class(x), " is an invalid input for entities")
}

init.VariableGroup <- function (.Object, group, entities, url.base=NULL, duplicates=FALSE, ...) {
    dots <- list(...)
    if ("variables" %in% names(dots)) entities <- dots$variables
    if ("name" %in% names(dots)) group <- dots$name
    .Object@group <- group
    .Object@entities <- .initEntities(entities, url.base)
    duplicates(.Object) <- duplicates
    return(.Object)
}
setMethod("initialize", "VariableGroup", init.VariableGroup)

##' @export
as.list.VariableOrder <- function (x, ...) x@graph

##' @export
as.list.VariableGroup <- function (x, ...) x@entities

##' Length of VariableOrder
##' @param x a VariableOrder
##' @return Integer: the number of VariableGroups in the VariableOrder
##' @name VariableOrder-length
NULL

##' @rdname VariableOrder-length
##' @export
setMethod("length", "VariableOrder", function (x) length(entities(x)))

##' Extract and update in VariableOrder and VariableGroup
##'
##' @param x a VariableOrder or VariableGroup
##' @param i an index. Numeric and logical indexing supported for both classes;
##' character indexing supported for VariableOrder, matching on VariableGroup
##' names
##' @param name Same as i but for \code{$}
##' @param j Invalid
##' @param value For update methods, an object equivalent in class to what is
##' being updated
##' @param ... additional arguments
##' @param drop Ignored
##' @return \code{[[} and \code{$} on a VariableOrder return the VariableGroup.
##' \code{[[} on VariableGroup returns the entity within, either a character
##' (URL) or nested VariableGroup. \code{[} and assignment methods return
##' objects of the same class as \code{x}
##' @name VariableOrder-extract
##' @aliases VariableOrder-extract
NULL

###############################
# 1. Extract from VariableOrder
###############################

##' @rdname VariableOrder-extract
##' @export
setMethod("[", c("VariableOrder", "ANY"), function (x, i, ..., drop=FALSE) {
    x@graph <- x@graph[i]
    return(x)
})
##' @rdname VariableOrder-extract
##' @export
setMethod("[", c("VariableOrder", "character"), function (x, i, ..., drop=FALSE) {
    w <- match(i, names(x))
    if (any(is.na(w))) {
        halt("Undefined groups selected: ", serialPaste(i[is.na(w)]))
    }
    callNextMethod(x, w, ..., drop=drop)
})

##' @rdname VariableOrder-extract
##' @export
setMethod("[[", c("VariableOrder", "ANY"), function (x, i, ...) {
    x@graph[[i]]
})

##' @rdname VariableOrder-extract
##' @export
setMethod("[[", c("VariableOrder", "character"), function (x, i, ...) {
    w <- match(i, names(x))
    callNextMethod(x, w, ..., drop=drop)
})

##' @rdname VariableOrder-extract
##' @export
setMethod("$", "VariableOrder", function (x, name) x[[name]])

###############################
# 2. Assign into VariableOrder
###############################

##' @rdname VariableOrder-extract
##' @export
setMethod("[<-", c("VariableOrder", "character", "missing", "VariableOrder"),
    function (x, i, j, value) {
        w <- match(i, names(x))
        if (any(is.na(w))) {
            halt("Undefined groups selected: ", serialPaste(i[is.na(w)]))
        }
        callNextMethod(x, w, value=value)
    })
##' @rdname VariableOrder-extract
##' @export
setMethod("[<-", c("VariableOrder", "ANY", "missing", "VariableOrder"),
   function (x, i, j, value) {
       x@graph[i] <- value@graph
       ## Ensure duplicates setting persists
       duplicates(x) <- duplicates(x)
       return(x)
   })


##' @rdname VariableOrder-extract
##' @export
setMethod("[[<-", c("VariableOrder", "character", "missing", "VariableGroup"),
    function (x, i, j, value) {
        w <- match(i, names(x))
        if (any(is.na(w))) {
            halt("Undefined group selected: ", serialPaste(i[is.na(w)]))
        }
        ## NextMethod: c("VariableOrder", "ANY", "missing", "VariableGroup")
        callNextMethod(x, w, value=value)
    })


.setNestedGroupByName <- function (x, i, j, value) {
    w <- match(i, names(x))
    value <- .initEntities(value)
    if (!duplicates(x)) {
        x <- setdiff_entities(x, value)
    }
    if (any(is.na(w))) {
        ## New group.
        entities(x) <- c(entities(x), VariableGroup(name=i, entities=value))
    } else {
        ## Existing group. Assign entities
        entities(x[[w]]) <- value
    }
    ## Ensure duplicates setting persists
    duplicates(x) <- duplicates(x)
    return(removeMissingEntities(x))
}
##' @rdname VariableOrder-extract
##' @export
setMethod("[[<-", c("VariableOrder", "character", "missing", "CrunchDataset"),
    .setNestedGroupByName)

##' @rdname VariableOrder-extract
##' @export
setMethod("[[<-", c("VariableOrder", "character", "missing", "VariableOrder"),
    function (x, i, j, value) {
        .setNestedGroupByName(x, i, j, entities(value))
    })

##' @rdname VariableOrder-extract
##' @export
setMethod("[[<-", c("VariableOrder", "ANY", "missing", "VariableGroup"),
    function (x, i, j, value) {
        if (!duplicates(x) && length(entities(value))) {
            x <- setdiff_entities(x, value)
        }
        x@graph[[i]] <- value
        ## Ensure duplicates setting persists
        duplicates(x) <- duplicates(x)
        return(removeMissingEntities(x))
    })

##' @rdname VariableOrder-extract
##' @export
setMethod("[[<-", c("VariableOrder", "ANY", "missing", "ANY"),
    function (x, i, j, value) {
        halt("Cannot assign an object of class ", dQuote(class(value)),
            " into a VariableOrder")
    })

##' @rdname VariableOrder-extract
##' @export
setMethod("[[<-", c("VariableOrder", "ANY", "missing", "NULL"),
    function (x, i, j, value) {
        x@graph[[i]] <- value
        return(x)
    })
##' @rdname VariableOrder-extract
##' @export
setMethod("[[<-", c("VariableOrder", "character", "missing", "NULL"),
    function (x, i, j, value) {
        w <- match(i, names(x))
        if (any(is.na(w))) {
            halt("Undefined group selected: ", serialPaste(i[is.na(w)]))
        }
        callNextMethod(x, w, value=value)
    })

##' @rdname VariableOrder-extract
##' @export
setMethod("$<-", "VariableOrder", function (x, name, value) {
    x[[name]] <- value
    return(x)
})

###############################
# 3. Extract from VariableGroup
###############################

##' @rdname VariableOrder-extract
##' @export
setMethod("[", c("VariableGroup", "ANY"), function (x, i, ..., drop=FALSE) {
    x@entities <- x@entities[i]
    return(x)
})
##' @rdname VariableOrder-extract
##' @export
setMethod("[", c("VariableGroup", "character"), function (x, i, ..., drop=FALSE) {
    w <- match(i, names(x))
    if (any(is.na(w))) {
        halt("Undefined groups selected: ", serialPaste(i[is.na(w)]))
    }
    callNextMethod(x, w, ..., drop=drop)
})


##' @rdname VariableOrder-extract
##' @export
setMethod("[[", c("VariableGroup", "character"), function (x, i, ...) {
    w <- match(i, names(x))
    if (any(is.na(w))) {
        halt("Undefined groups selected: ", serialPaste(i[is.na(w)]))
    }
    callNextMethod(x, w, ..., drop=drop)
})

##' @rdname VariableOrder-extract
##' @export
setMethod("[[", c("VariableGroup", "ANY"), function (x, i, ...) {
    x@entities[[i]]
})

##' @rdname VariableOrder-extract
##' @export
setMethod("$", "VariableGroup", function (x, name) x[[name]])

###############################
# 4. Assign into VariableGroup
###############################

##' @rdname VariableOrder-extract
##' @export
setMethod("[[<-", c("VariableGroup", "character", "missing", "CrunchDataset"),
    .setNestedGroupByName)

##' @rdname VariableOrder-extract
##' @export
setMethod("[[<-", c("VariableGroup", "character", "missing", "list"),
    .setNestedGroupByName)

##' @rdname VariableOrder-extract
##' @export
setMethod("[[<-", c("VariableGroup", "character", "missing", "character"),
    .setNestedGroupByName)

##' @rdname VariableOrder-extract
##' @export
setMethod("[[<-", c("VariableGroup", "character", "missing", "VariableOrder"),
    function (x, i, j, value) {
        .setNestedGroupByName(x, i, j, entities(value))
    })

##' @rdname VariableOrder-extract
##' @export
setMethod("[[<-", c("VariableGroup", "character", "missing", "VariableGroup"),
    function (x, i, j, value) {
        .setNestedGroupByName(x, i, j, entities(value))
    })

##' @rdname VariableOrder-extract
##' @export
setMethod("[[<-", c("VariableGroup", "ANY", "missing", "VariableGroup"),
    function (x, i, j, value) {
        entities(x)[[i]] <- value
        return(x)
    })

##' @rdname VariableOrder-extract
##' @export
setMethod("$<-", "VariableGroup", function (x, name, value) {
    x[[name]] <- value
    return(x)
})

##' Get un(grouped) VariableGroups
##'
##' "ungrouped" is a magic VariableGroup that contains all variables not found
##' in groups at a given level of nesting.
##' @param var.order an object of class VariableOrder or VariableGroup
##' @return For grouped(), a VariableOrder/Group, respectively, with "ungrouped"
##' omitted. For ungrouped(), a VariableGroup.
##' @seealso \code{\link{VariableOrder}}
##' @export
grouped <- function (var.order) {
    Filter(Negate(is.character), var.order)
}

##' @rdname grouped
##' @export
ungrouped <- function (var.order) {
    return(VariableGroup(name="ungrouped",
        entities=entities(Filter(is.character, var.order))))
}

setdiff_entities <- function (x, ents, remove.na=FALSE) {
    ## Remove "ents" (variable references) anywhere they appear in x (Order)
    if (!is.character(ents)) {
        ## Get just the entity URLs
        ents <- entities(ents, simplify=TRUE)
    }

    if (inherits(x, "VariableOrder") || inherits(x, "VariableGroup")) {
        entities(x) <- setdiff_entities(entities(x), ents)
    } else if (is.list(x)) {
        ## We're inside entities, which may have nested groups
        grps <- vapply(x, inherits, logical(1), what="VariableGroup")
        x[grps] <- lapply(x[grps], setdiff_entities, ents)
        matches <- unlist(x[!grps]) %in% ents
        if (any(matches)) {
            ## Put in NAs so that any subsequent assignment into this object
            ## assigns into the right position. Then strip NAs after
            x[!grps][matches] <- rep(list(NA_character_), sum(matches))
        }
    }
    if (remove.na) {
        x <- removeMissingEntities(x)
    }
    return(x)
}

removeMissingEntities <- function (x) {
    ## Remove NA entries, left by setdiff_entities, from @graph/entities
    if (inherits(x, "VariableOrder") || inherits(x, "VariableGroup")) {
        entities(x) <- removeMissingEntities(entities(x))
    } else if (is.list(x)) {
        ## We're inside entities, which may have nested groups
        grps <- vapply(x, inherits, logical(1), what="VariableGroup")
        x[grps] <- lapply(x[grps], removeMissingEntities)
        drops <- vapply(x[!grps], is.na, logical(1))
        if (any(drops)) {
            x <- x[-which(!grps)[drops]]
        }
    }
    return(x)
}
