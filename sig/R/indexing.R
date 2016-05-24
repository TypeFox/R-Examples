#' Indexing for siglists
#' 
#' Get or set a subset of a siglist.
#' @param x A \code{siglist} object.
#' @param i An integer vector index.
#' @param ... Passed from other index methods.
#' @param value A value to set the subset to.
#' @return A siglist.
#' @seealso \code{\link[base]{Extract}}
#' @examples
#' methods_sigs <- list_sigs(pkg2env(methods))
#' methods_sigs[1:5]
#' methods_sigs[[1]]
#' @name [
#' @aliases Extract.siglist
#' @aliases [.siglist
#' @aliases [[.siglist
#' @aliases [<-.siglist
#' @aliases [[<-.siglist
#' @rdname Extract.siglist
#' @method [ siglist
#' @export
`[.siglist` <- function(x, i, ...)
{
  structure(
    unclass(x)[i], 
    class = "siglist"
  )
}

#' @rdname Extract.siglist
#' @method [[ siglist
#' @export
`[[.siglist` <- function(x, i, ...)
{
  structure(
    unclass(x)[[i]], 
    class = "sig"
  )
}

#' @rdname Extract.siglist
#' @method [<- siglist
#' @export
`[<-.siglist` <- function(x, ..., value)
{
  x <- unclass(x)
  value <- lapply(value, as.sig)
  x[...] <- value
  structure(
    x, 
    class = "siglist"
  )
}

#' @rdname Extract.siglist
#' @method [[<- siglist
#' @export
`[[<-.siglist` <- function(x, ..., value)
{
  x <- unclass(x)
  value <- as.sig(value)
  x[...] <- value
  structure(
    x, 
    class = "siglist"
  )
}
