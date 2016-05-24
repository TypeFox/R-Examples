##' Coerce an object to a mutaframe. Supported types include
##' \code{data.frame}, or anything coercible to one.
##'
##' @title Coercion to mutaframe
##' @param x the object to coerce
##' @param ... arguments passed to methods
##' @return a mutaframe
##' @rdname as.mutaframe
##' @export
as.mutaframe <- function(x, ...) UseMethod("as.mutaframe")

##' @rdname as.mutaframe
##' @method as.mutaframe mutaframe
##' @export
as.mutaframe.mutaframe <- function(x, ...) {
  cl <- oldClass(x)
  i <- match("mutaframe", cl)
  if (i > 1L) 
    class(x) <- cl[-(1L:(i - 1L))]
  x
}

##' @rdname as.mutaframe
##' @method as.mutaframe data.frame
##' @export
as.mutaframe.data.frame <- function(x, ...) .mutaframe(x, rownames(x))

##' @rdname as.mutaframe
##' @method as.mutaframe default
##' @export
as.mutaframe.default <- function(x, ...) as.mutaframe(as.data.frame(x, ...))

##' Coerces a mutaframe to a \code{data.frame}
##'
##' @title Coercion to data.frame
##' @param x a mutaframe
##' @param row.names character vector of rownames, defaults to
##' rownames of \code{x}
##' @param optional see \code{\link{as.data.frame}}
##' @param ... see \code{\link{as.data.frame}}
##' @return a \code{data.frame}
##' @method as.data.frame mutaframe
##' @export
as.data.frame.mutaframe <- function(x, row.names = rownames(x),
                                    optional = FALSE, ...)
{
  cols <- lapply(names(x), function(j) x[[j]])
  names(cols) <- names(x)
  df <- as.data.frame(cols, optional = optional, ...)
  ## we set row.names this way for cases where we have no columns, but >0 rows
  ## as.data.frame complains in that case
  attr(df, "row.names") <- row.names
  df
}

##' Tests whether an object is a \code{mutaframe}
##'
##' @title Test for mutaframes
##' @param x an object to check
##' @return \code{TRUE} if \code{x} is an instance of a class that
##' inherits from \code{mutaframe}; otherwise, \code{FALSE}
##' @export
is.mutaframe <- function(x) inherits(x, "mutaframe")

##' Coerces a mutaframe to a list
##'
##' @title Coercion to list
##' @param x a mutaframe
##' @param ... ignored
##' @return a list, with one element for each mutaframe column
##' @method as.list mutaframe
##' @export
as.list.mutaframe <- function(x, ...) lapply(x, do.call, list())
## as.list.environment does not resolve active bindings
