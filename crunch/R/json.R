##' @rdname tojson-crunch
##' @export
setMethod("jsonprep", "Categories", function (x, ...) jsonprep(I(x@.Data)))

##' @rdname tojson-crunch
##' @export
setMethod("jsonprep", "list", function (x, ...) lapply(x, jsonprep, ...))

##' @rdname tojson-crunch
##' @export
setMethod("jsonprep", "ANY", function (x, ...) x)


.jsonprep.vargroup <- function (x, ...) {
    ents <- x@entities
    if (length(ents) == 0) {
        ## toJSON(character(0)) is [""], which is length 1 :(
        ents <- list() ## but toJSON(list()) is []
    } else if (is.list(ents)) {
        nested.groups <- vapply(ents, inherits, logical(1),
            what="VariableGroup")
        ents[nested.groups] <- lapply(ents[nested.groups], .jsonprep.vargroup)
    }
    return(structure(list(I(ents)), .Names=x@group))
}

##' @rdname tojson-crunch
##' @export
setMethod("jsonprep", "VariableOrder",
    function (x, ...) jsonprep(list(graph=x@graph, ...)))

##' @rdname tojson-crunch
##' @export
setMethod("jsonprep", "VariableGroup", .jsonprep.vargroup)


##' @importFrom jsonlite toJSON
##' @rdname tojson-crunch
##' @export
toJSON <- function (x, ...) {
    out <- jsonlite::toJSON(jsonprep(x), auto_unbox=TRUE, null="null", na="null", force=TRUE, ...)
    # cat(out)
    return(out)
}
