#--------------------------------------------------------------------
#   grid.par.R (npsp package)
#--------------------------------------------------------------------
#   grid.par    S3 class and methods
#       as.data.frame.grid.par(x, row.names = names(x$min), optional = FALSE)
#       print.grid.par(x)
#       dimnames.grid.par(x)
#       dim.grid.par(x)
#   coords      S3 generic
#       coords.grid.par(x)
#   coordvalues S3 generic
#       coordvalues.grid.par(x)
#
# PENDENTE:
#   - renombrar coords -> coordinates ? (problemas có método S4 de sp?)
#   - is.grid.par()
#   - ndim()
#
#   (c) R. Fernandez-Casal         Last revision: Aug 2012
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# grid.par(n, min, max = min + (n-1)*lag, lag = (max-min)/(n-1)) 
# 'Based' on GridTopology-class of sp package
#--------------------------------------------------------------------
#' Grid parameters (S3 class "grid.par")
#'
#' Defines a full regular (spatial) grid. 
#' Constructor function of the \code{grid.par}-\code{\link{class}}.
#' @aliases grid.par-class
#' @param  n integer vector; number of nodes in each dimension. 
#' @param  min vector; minimun coordinates values.
#' @param  max vector; maximum coordinates values (optional).
#' @param  lag vector; lag in each dimension (optional). 
#' @param  dimnames character vector; names used to label the dimensions.
#' @details All parameters must have the same length. 
#' Only one of the arguments \code{max} or \code{lag} must be specified. 
#' 
#' S3 'version' of the \code{\link[sp]{GridTopology}}-\code{\link{class}} of the \pkg{sp} package.
#' @return Returns an object of class \code{grid.par}, a list with the arguments as components 
#'    and an additional component \code{$nd = length(n)}.
#' @seealso \code{\link{data.grid}}.
#' @export
#' @examples
#' grid.par(n = c(100, 100), min = c(-10, 42), max = c(-7.5, 44))
#' grid.par(n = c(100, 100), min = c(-10, 42), lag = c(0.03, 0.02))
grid.par <- function(n, min, max = min + (n-1)*lag, lag = (max-min)/(n-1), dimnames = names(min)) {
#
# PENDENTE:
#   - pasar a metodo S3 y esta a default
#   - ns <- names(x$min)
#--------------------------------------------------------------------
    if ( !missing(lag) && !missing(max) ) {
        warning("argument 'lag' ignored (set to default value)")
        lag <- (max-min)/(n-1)
    }     
    nd <- length(n)
    if (length(min) != nd) stop("argument 'min' has incorrect dimension")
    if (length(max) != nd) stop("argument 'max' has incorrect dimension")
    if (length(lag) != nd) stop("argument 'lag' has incorrect dimension")
    if (is.null(dimnames)) dimnames <- paste("x", seq_len(nd), sep="")
    result <- list(nd = nd, n = n, min = min, max = max, lag = lag, dimnames = dimnames)
    oldClass(result) <- "grid.par"
    return(result)
#--------------------------------------------------------------------
} # grid.par






