#--------------------------------------------------------------------
#   npsp-internals.R (npsp package)
#--------------------------------------------------------------------
#   grid.par  S3 methods
#       print.grid.par(x, ...)
#       dim.grid.par(x)
#       dimnames.grid.par(x)
#       as.data.frame.grid.par(x, row.names, optional, ...)
#   data.grid S3 methods
#       dim.data.grid(x)
#       dimnames.data.grid(x)      
#   is.data.grid(x)       
#
#   (c) R. Fernandez-Casal         Last revision: Aug 2013
#--------------------------------------------------------------------



#--------------------------------------------# grid.par #--------------------------------------------####


#--------------------------------------------------------------------
#' @name npsp-internals
#' @aliases print.grid.par
#' @title npsp internal and secondary functions
#' @description Listed below are supporting functions for the major methods in npsp.
#' @method print grid.par
#' @keywords internal
#' @export
print.grid.par <- function(x, ...) {
#--------------------------------------------------------------------
  result <- as.data.frame(x)
	print(result)
	invisible(result)
}

#--------------------------------------------------------------------
#' @rdname npsp-internals
#' @method dim grid.par
#' @keywords internal
#' @export
dim.grid.par <- function(x) x$n
#--------------------------------------------------------------------

#--------------------------------------------------------------------
#' @rdname npsp-internals
#' @method dimnames grid.par
#' @keywords internal
#' @export
dimnames.grid.par <- function(x) {
#--------------------------------------------------------------------
    if (!inherits(x, "grid.par"))
      stop("function only works for objects of class (or extending) 'grid.par'")
    ns <- x$dimnames
    if (is.null(ns)) ns <- paste("x",seq_len(x$nd),sep="")
    return(ns)
}

#--------------------------------------------------------------------
# as.data.frame.grid.par(x, row.names = names(x$min), optional = FALSE, ...){
#' @rdname npsp-internals
#' @method as.data.frame grid.par
#' @keywords internal
#' @export
as.data.frame.grid.par <- function(x, row.names = dimnames(x), optional = FALSE, ...){
#--------------------------------------------------------------------
    if (is.null(row.names)) row.names <- paste("x",seq_len(x$nd),sep="")
    data.frame(n = x$n, min = x$min, max = x$max, lag = x$lag, row.names = row.names)
}


#--------------------------------------------# data.grid #--------------------------------------------####



#--------------------------------------------------------------------
#' @rdname npsp-internals
#' @keywords internal
#' @export
is.data.grid <- function(x) inherits(x, "data.grid")
#--------------------------------------------------------------------


#--------------------------------------------------------------------
#' @rdname npsp-internals
#' @method dim data.grid
#' @keywords internal
#' @export
dim.data.grid <- function(x) dim.grid.par(x$grid)
#--------------------------------------------------------------------


#--------------------------------------------------------------------
#' @rdname npsp-internals
#' @method dimnames data.grid
#' @keywords internal
#' @export
dimnames.data.grid <- function(x) dimnames.grid.par(x$grid)
#--------------------------------------------------------------------




