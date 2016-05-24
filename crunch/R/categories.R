setValidity("Categories", function (object) {
    are.cats <- vapply(object, is.category, logical(1))
    if (!all(are.cats)) {
        badcount <- sum(!are.cats)
        return(paste0("Invalid categories: ", badcount,
            ifelse(badcount>1,
                " elements are not Crunch category objects.",
                " element is not a Crunch category object.")))
    }
    if (any(duplicated(names(object)))) {
        return("Invalid category names: must be unique")
    }
    if (any(duplicated(ids(object)))) {
        return("Invalid category ids: must be unique")
    }
    return(TRUE)
})

setMethod("initialize", "Categories", function (.Object, ...) {
    .Object@.Data <- lapply(..1,
        function (x) try(Category(data=x), silent=TRUE))
    validObject(.Object)
    return(.Object)
})

is.categories <- function (x) inherits(x, "Categories")

concatenateCategories <- function (...) {
    ## c() S3 method for categories. Dispatch is on ..1
    dots <- list(...)
    iscat <- vapply(dots, is.category, logical(1))
    iscats <- vapply(dots, is.categories, logical(1))
    if (!all(iscat | iscats)) {
        stop("Invalid categories")
    }
    dots[iscat] <- lapply(dots[iscat], function (x) list(x))
    dots[iscats] <- lapply(dots[iscats], function (x) x@.Data)
    return(Categories(data=do.call(c, dots)))
}

##' S3 method to concatenate Categories and Category objects
##'
##' @param ... see \code{\link[base]{c}}
##' @return An object of class \code{\link{Categories}}
##' @name c-categories
##' @export
##' @examples
##' cat.a <- Category(name="First", id=1, numeric_value=1, missing=FALSE)
##' cat.b <- Category(name="Second", id=2)
##' cat.c <- Category(name="Third", id=3, missing=TRUE)
##' cats.1 <- Categories(cat.a, cat.b)
##' identical(cats.1, c(cat.a, cat.b))
##' identical(c(cats.1, cat.c), Categories(cat.a, cat.b, cat.c))
c.Categories <- concatenateCategories

##' @rdname c-categories
##' @export
c.Category <- concatenateCategories

##' @rdname Categories
##' @export
setMethod("[", c("Categories", "ANY"), function (x, i, ...) {
    x@.Data <- x@.Data[i]
    return(x)
})

##' @rdname Categories
##' @export
setMethod("[", c("Categories", "numeric"), function (x, i, ...) {
    invalid.indices <- setdiff(abs(i), seq_along(x@.Data))
    if (length(invalid.indices)) {
        halt("subscript out of bounds: ", serialPaste(invalid.indices))
    }
    x@.Data <- x@.Data[i]
    return(x)
})

##' @rdname Categories
##' @export
setMethod("[<-", c("Categories", "ANY"), function (x, i, ..., value) {
    x@.Data[i] <- Categories(data=value)
    return(x)
})

##' @rdname Categories
##' @export
setMethod("names", "Categories", function (x) vapply(x, name, character(1)))

##' @rdname Categories
##' @export
setMethod("values", "Categories", function (x) vapply(x, value, numeric(1)))

##' @rdname Categories
##' @export
setMethod("ids", "Categories", function (x) vapply(x, id, integer(1)))

setNames <- function (x, value) {
    if (is.null(value) || !is.character(value)) {
        halt('Names must be of class "character"')
    }
    if (!identical(length(x), length(value))) {
        halt("Invalid names: supplied ", length(value), " names for ", length(x), " categories")
    }
    x[] <- mapply(setName, x, value=value, SIMPLIFY=FALSE)
    return(x)
}
setValues <- function (x, value) {
    x[] <- mapply(setValue, x[], value=value, SIMPLIFY=FALSE)
    return(x)
}

##' @rdname Categories
##' @export
setMethod("names<-", "Categories", setNames)
##' @rdname Categories
##' @export
setMethod("values<-", "Categories", setValues)
##' @rdname Categories
##' @export
setMethod("ids<-", "Categories", function (x, value) {
    if (!identical(ids(x), value)) {
        halt("Cannot modify category ids")
    }
    return(x)
})

.na.omit.categories <- function (object, ...) {
    missings <- vapply(object, function (x) isTRUE(x$missing), logical(1),
        USE.NAMES=FALSE)
    if (any(missings)) {
        object <- object[!missings]
        attr(object, "na.action") <- which(missings)
        attr(object, "class") <- "omit"
    }
    return(object)
}

##' Omit missing categories
##' @param object Categories
##' @param ... additional arguments, ignored
##' @return \code{object} with any categories that have missing: TRUE excluded
##' @name na-omit-categories
NULL

##' @rdname na-omit-categories
##' @export
setMethod("na.omit", "Categories", function (object, ...) {
    Categories(data=.na.omit.categories(object))
})

##' is.na for Categories
##'
##' @param x Categories or a single Category
##' @param value To change the missingness of categories, supply either (1)
##' a logical vector of equal length of the categories (or length 1 for the
##' Category method), or (2) the names of the categories to mark as missing.
##' If supplying the latter, any categories already indicated as missing will
##' remain missing.
##' @return Getters return logical, a named vector in the case of the Categories
##' method; setters return \code{x} duly modified.
##' @name is-na-categories
NULL

##' @rdname is-na-categories
##' @aliases is-na-categories
##' @export
setMethod("is.na", "Categories", function (x) structure(vapply(x, is.na, logical(1), USE.NAMES=FALSE), .Names=names(x)))

n2i <- function (x, cats, strict=TRUE) {
    ## Convert x from category names to the corresponding category ids
    out <- ids(cats)[match(x, names(cats))]
    if (strict && any(is.na(out))) {
        halt(ifelse(sum(is.na(out)) > 1, "Categories", "Category"),
            " not found: ", serialPaste(dQuote(x[is.na(out)])))
    }
    return(out)
}

i2n <- function (x, cats, strict=TRUE) {
    ## Convert x from category ids to the corresponding category names
    out <- names(cats)[match(x, ids(cats))]
    if (strict && any(is.na(out))) {
        halt(ifelse(sum(is.na(out)) > 1, "Categories", "Category"),
            " not found: ", serialPaste(dQuote(x[is.na(out)])))
    }
    return(out)
}

##' @rdname is-na-categories
##' @export
setMethod("is.na<-", c("Categories", "character"), function (x, value) {
    ix <- match(value, names(x))
    if (any(is.na(ix))) {
        halt(ifelse(sum(is.na(ix)) > 1, "Categories", "Category"),
            " not found: ", serialPaste(dQuote(value[is.na(ix)])))
    }
    x[ix] <- lapply(x[ix], `is.na<-`, value=TRUE)
    return(x)
})

##' @rdname is-na-categories
##' @export
setMethod("is.na<-", c("Categories", "logical"), function (x, value) {
    stopifnot(length(x) == length(value))
    x@.Data <- mapply(function (x, value) {
            is.na(x) <- value
            return(x)
        }, x=x@.Data, value=value, USE.NAMES=FALSE, SIMPLIFY=FALSE)
    return(x)
})

setMethod("lapply", "Categories", function (X, FUN, ...) {
    X@.Data <- lapply(X@.Data, FUN, ...)
    return(X)
})
