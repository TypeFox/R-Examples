##' Indicate how categories represent a dichotomized value
##'
##' Multiple Response variables are essentially Categorical Arrays that have
##' had a category or categories indicated as the "selected" value. These
##' methods let you set that state.
##'
##' \code{dichotomize} lets you specify which categories are "selected", while
##' \code{undichotomize} strips that selection information. Dichotomize converts
##' a Categorical Array to a Multiple Response, and undichotomize converts back.
##'
##' @param x Categories or a Variable subclass that has Categories
##' @param i For the \code{dichotomize} methods, the numeric or logical indices
##' of the categories to mark as "selected", or if character, the Category
##' "names". Note that unlike some other categorical variable methods,
##' numeric indices are positional, not with reference to category ids.
##' @return Categories or the Variable, (un)dichotomized accoringly
##' @name dichotomize
##' @aliases dichotomize is.dichotomized undichotomize
##' @seealso \code{\link{describe-category}}
NULL

##' @rdname dichotomize
##' @export
setMethod("is.dichotomized", "Categories",
    function (x) any(vapply(x, is.selected, logical(1))))

.dichotomize.categories <- function (x, i) {
    ## Internal method for dichtomizing Categories (or lists)
    x[i] <- lapply(x[i], function (a) {
        a$selected <- TRUE
        return(a)
    })
    return(x)
}

##' @rdname dichotomize
##' @export
setMethod("dichotomize", c("Categories", "numeric"), .dichotomize.categories)
##' @rdname dichotomize
##' @export
setMethod("dichotomize", c("Categories", "logical"), .dichotomize.categories)
##' @rdname dichotomize
##' @export
setMethod("dichotomize", c("Categories", "character"), function (x, i) {
    ind <- names(x) %in% i
    if (!any(ind)) {
        halt("Category not found") ## make nicer error message
    }
    return(dichotomize(x, ind))
})

##' @rdname dichotomize
##' @export
setMethod("undichotomize", "Categories", function (x) {
    x[] <- lapply(x[], function (a) {
        a$selected <- FALSE
        return(a)
    })
    return(x)
})

.dichotomize.var <- function (x, i) {
    newcats <- dichotomize(categories(x), i)
    categories(x) <- newcats
    if (is.dichotomized(newcats)) {
        ## Do this to avoid needing to refresh the variable catalog
        x@tuple@body$type <- "multiple_response"
    }
    invisible(CrunchVariable(tuple(x)))
}
.undichotomize.var <- function (x) {
    categories(x) <- undichotomize(categories(x))
    ## Do this to avoid needing to refresh the variable catalog
    x@tuple@body$type <- "categorical_array"
    invisible(CrunchVariable(tuple(x)))
}

##' @rdname dichotomize
##' @export
setMethod("dichotomize", "CategoricalVariable", .dichotomize.var)
##' @rdname dichotomize
##' @export
setMethod("dichotomize", "CategoricalArrayVariable", .dichotomize.var)
##' @rdname dichotomize
##' @export
setMethod("undichotomize", "CategoricalVariable", .undichotomize.var)
##' @rdname dichotomize
##' @export
setMethod("undichotomize", "CategoricalArrayVariable", .undichotomize.var)
