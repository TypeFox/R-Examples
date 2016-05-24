#' @export
dimnames.spam <- function(x) attr(x, 'spamdimnames')
#' @export
`dimnames<-.spam` <- function(x, value) {attr(x, 'spamdimnames') <- value; x}

#' Extract parts of a sparse \code{spam} matrix
#'
#' @name [
#' @param x object to extract from.
#' @param i row identifiers.
#' @param j column identifiers.
#' @param drop logical indicating that dimensions should be dropped.
#' @param ... additional arguments.
#' @aliases [,spam,character,character,logical-method
#' @docType methods
#' @rdname extract-methods
#' @export
setMethod("[", signature(x = "spam", i="character", j="character", drop = "logical"), function(x, i, j, ..., drop) return(GetElements.spam(x, i, j, drop)))

#' @name [
#' @aliases [,spam,character,character,missing-method
#' @docType methods
#' @rdname extract-methods
#' @export
setMethod("[", signature(x = "spam", i="character", j="character", drop = "missing"), function(x, i, j, ..., drop) return(GetElements.spam(x, i, j, drop)))

#' @name [
#' @aliases [,spam,character,missing,logical-method
#' @docType methods
#' @rdname extract-methods
#' @export
setMethod("[", signature(x = "spam", i="character", j="missing", drop = "logical"), function(x, i, j, ..., drop) return(GetElements.spam(x, i, j, drop)))

#' @name [
#' @aliases [,spam,character,missing,missing-method
#' @docType methods
#' @rdname extract-methods
#' @export
setMethod("[", signature(x = "spam", i="character", j="missing", drop = "missing"), function(x, i, j, ..., drop) return(GetElements.spam(x, i, j, drop)))

#' @name [
#' @aliases [,spam,missing,character,logical-method
#' @docType methods
#' @rdname extract-methods
#' @export
setMethod("[", signature(x = "spam", i="missing", j="character", drop = "logical"), function(x, i, j, ..., drop) return(GetElements.spam(x, i, j, drop)))

#' @name [
#' @aliases [,spam,missing,character,missing-method
#' @docType methods
#' @rdname extract-methods
#' @export
setMethod("[", signature(x = "spam", i="missing", j="character", drop = "missing"), function(x, i, j, ..., drop) return(GetElements.spam(x, i, j, drop)))


GetElements.spam <- function(x, i, j, drop) {
  if(!missing(i)) {
    if (is.character(i)){
      if (is.null(rownames(x))) stop("Row names do not exist.")
      else i <- match(i, rownames(x))
    } else stop("I do not understand the subset criteria!")
  }

  if(!missing(j)) {
    if (is.character(j)){
      if (is.null(colnames(x))) stop("Column names do not exist.")
      else j <- match(j, colnames(x))
    } else stop("I do not understand the subset criteria!")
  }

  if(missing(i)) x[seq_len(nrow(x)), j, drop=ifelse(missing(drop), FALSE, drop)]
  else if(missing(j)) x[i, , drop=ifelse(missing(drop), TRUE, drop)]
  else x[i, j, drop=ifelse(missing(drop), TRUE, drop)]
}
