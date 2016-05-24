is.category <- function (x) inherits(x, "Category")

setValidity("Category", function (object) {
    is.cat <- all(c("id", "name") %in% names(object))
    if (!all(is.cat)) {
        val <- "Not a category"
    } else {
        val <- TRUE
    }
    return(val)
})

init.Category <- function (.Object, ...) {
    .Object <- callNextMethod()
    ## Make sure category elements are sorted so that identical categories are
    ## evaluated identically. Order doesn't matter for object, but R lists are
    ## ordered.
    s <- order(.Object@names)
    .Object@.Data <- .Object@.Data[s]
    .Object@names <- .Object@names[s]
    return(.Object)
}
setMethod("initialize", "Category", init.Category)

.no.data <- list(
    id=-1L,
    name="No Data",
    numeric_value=NULL,
    missing=TRUE
)

setName <- function (x, value) {
    x[["name"]] <- value
    return(x)
}
setValue <- function (x, value) {
    value_to_set <- suppressWarnings(as.numeric(value))
    if (is.na(value_to_set) && !is.na(value)) {
        halt("Category values must be numeric")
    }
    x[["numeric_value"]] <- value_to_set
    return(x)
}

##' Access Category fields directly
##'
##' Don't do this. Instead, use the category setters.
##'
##' @param x a Category
##' @param name a field within \code{x}
##' @param value a value for that field to update
##' @return \code{$} returns the value of the desired field. Setter
##' returns \code{x} duly modified.
##' @seealso \code{\link{describe-category}}
##' @name category-extract
NULL

##' @rdname category-extract
##' @export
setMethod("$", "Category", function (x, name) x[[name]])
##' @rdname category-extract
##' @export
setMethod("$<-", "Category", function (x, name, value) {
    x[[name]] <- value
    return(x)
})

##' Category attributes
##'
##' @param x a Category
##' @param value For the setters, an appropriate value to set
##' @return \code{name} returns character; \code{value} and \code{id} return
##' numeric; value but not id may be \code{NA}; \code{is.selected} returns
##' logical indicating whether this Category is a "selected" dichotomy. Setters
##' return \code{x} duly modified.
##' @name describe-category
##' @aliases value value<- id is.selected describe-category
##' @seealso \code{\link{Categories}} \code{\link{dichotomize}}
NULL

##' @rdname describe-category
##' @export
setMethod("name", "Category", function (x) x[["name"]])
##' @rdname describe-category
##' @export
setMethod("name<-", c("Category", "character"), setName)
##' @rdname describe-category
##' @export
setMethod("name<-", c("ANY", "ANY"), function (x, value) {
    halt('Names must be of class "character"')
})
##' @rdname describe-category
##' @export
setMethod("value", "Category", function (x) {
    v <- x[["numeric_value"]]
    return(ifelse(is.null(v), NA_real_, as.numeric(v)))
})
##' @rdname describe-category
##' @export
setMethod("value<-", "Category", setValue)
##' @rdname describe-category
##' @export
setMethod("id", "Category", function (x) as.integer(x[["id"]]))

##' @rdname describe-category
##' @export
setMethod("is.selected", "Category", function (x) isTRUE(x$selected))

##' @rdname is-na-categories
##' @export
setMethod("is.na", "Category", function (x) isTRUE(x$missing))

##' @rdname is-na-categories
##' @export
setMethod("is.na<-", c("Category", "logical"), function (x, value) {
    stopifnot(length(value) == 1)
    x$missing <- isTRUE(value)
    return(x)
})
