cubeDims <- function (cube) {
    ## Grab the row/col/etc. labels from the cube
    elements <- lapply(cube$result$dimensions, function (a) {
        ## If enumerated, will be "elements", not "categories"
        a$type$categories %||% a$type$elements
    })
    dimnames <- lapply(elements, function (d) {
        list(
            name=vapply(d, elementName, character(1)),
            any.or.none=vapply(d, elementIsAnyOrNone, logical(1)),
            missing=vapply(d, function (el) isTRUE(el$missing), logical(1))
        )
    })
    names(dimnames) <- vapply(cube$result$dimensions,
        function (a) a$references$alias, character(1))
    return(CubeDims(dimnames))
}

elementName <- function (el) {
    ## Given a category or element in a cube dim, what's its name?
    out <- el$value
    if (is.null(out)) {
        ## This is probably categorical. Try "name" instead of "value".
        out <- el$name
    } else if (is.list(out)) {
        if (length(out) == 2 && is.null(names(out))) {
            ## el$value is bin boundaries, as in a binned numeric.
            out <- paste(unlist(out), collapse="-")
        } else {
            ## This is probably a subvariable. Look for its name.
            out <- out$references$name
        }
    }
    if (is.null(out)) {
        ## Damn. You may be here because you're hitting missing values in an
        ## array or multiple response, or the __any__ or __none__ values.
        ## Bail out.
        out <- "<NA>"
    }
    out <- as.character(out)
    return(out)
}

elementIsAnyOrNone <- function (el) {
    is.list(el$value) && ## Element has $value and value is a list
        "id" %in% names(el$value) && ## "value" has names (is not bin)
        el$value$id %in% c("__any__", "__none__")
}

##' Methods on Cube objects
##'
##' These methods provide an \code{array}-like interface to the CrunchCube
##' object.
##'
##' @param x a CrunchCube or its CubeDims component.
##'
##' @return Generally, the same shape of result that each of these functions
##' return when applied to an \code{array} object.
##' @name cube-methods
##' @aliases cube-methods
##' @seealso \code{\link{cube-computing}}
NULL

##' @rdname cube-methods
##' @export
setMethod("dimnames", "CubeDims", function (x) {
    lapply(x, function (a) a$name)
})

##' @rdname cube-methods
##' @export
setMethod("dim", "CubeDims",
    function (x) vapply(dimnames(x), length, integer(1), USE.NAMES=FALSE))

##' @rdname cube-methods
##' @export
setMethod("is.na", "CubeDims", function (x) lapply(x, function (a) a$missing))

anyOrNone <- function (x) {
    lapply(x, function (a) a$any.or.none)
}
