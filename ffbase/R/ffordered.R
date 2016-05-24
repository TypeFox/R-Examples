#' Add the order of a ff vector of the \code{\link{ff}} vector \code{x}
#' 
#' Add a ff vector that contains the order of the \code{\link{ff}} vector \code{x} as an attribute.
#' The order can be retrieved using \code{ffordered}. Note that you have to assign the result to the
#' original vector \code{x}.
#' @example ../examples/ffordered.R
#' @param x \code{ff} vector to be indexed
#' @param addsorted should the sorted values also be stored in \code{ffsorted}?
#' @param ... parameters that will be passed on to \code{\link{fforder}}. 
#' @return The updated vector \code{x}
#' @rdname ffordered
#' @export
#' @aliases addfforder ffordered ffsorted
addfforder <- function(x, addsorted=FALSE, ...){
  ordered <- fforder(x)
  filename(ordered) <- sub("(\\.ff)?$", ".fforder", filename(x))
  attr(x, "fforder") <- ordered
  if (addsorted){
     sorted <- x[ordered]
     is.sorted(sorted) <- TRUE
     filename(sorted) <- sub("(\\.ff)?$", ".ffsort", filename(x))
     attr(x, "ffsort") <- sorted
  }
  x
}

#' @rdname ffordered
#' @export
ffordered <- function(x){
  ordered <- attr(x, "ffordered")
  if (is.null(ordered)){
    ordered <- fforder(x)
  } 
  ordered
}

#' @rdname ffordered
#' @export
ffsorted <- function(x){
  sorted <- attr(x, "ffsort")
  if (is.null(sorted)){
    sorted <- x[ffordered(x)]
    is.sorted(sorted) <- TRUE
  }
  sorted
}

# ### Quick testing
# x <- ff(runif(10))
# x <- addfforder(x)
# x[]
# ffordered(x)[]
# ffsorted(x)[]
