#--------------------------------------------------------------------
#   data.grid.R (npsp package)
#--------------------------------------------------------------------
#   data.grid   S3 class and methods
#       coords.data.grid(x)
#       coordvalues.data.grid(x)
#       dimnames.data.grid(x)
#       dim.data.grid(x)
#   interp.data.grid()
#
# PENDENTE:
#   - exemplos
#   - as.data.grid()
#   - as.data.frame.data.grid()
#       !sapply(data.grid, is.list) & (sapply(data.grid, length) == prod(dim(data.grid)))
#       coords = TRUE, ns <- names(coords) if(any(ns %in% names) ns <- paste("coord", ns, sep=".")
#
#   (c) R. Fernandez-Casal         Last revision: Aug 2012
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# data.grid(..., grid = NULL)
# 'Equivalent' to SpatialGridDataFrame-class of sp package
#--------------------------------------------------------------------
#' Gridded data (S3 class "data.grid")
#' 
#' Defines data on a full regular (spatial) grid. 
#' Constructor function of the \code{data.grid}-\code{\link{class}}.
#' @aliases data.grid-class
#' @param  ... vectors or arrays of data with length equal to \code{prod(grid$n)}. 
#' @param  grid a \code{\link{grid.par}}-\code{\link{class}} object (optional).
#' @details If parameter \code{grid.par} is not specified it is set from first argument.
#'
#' S3 "version" of the \code{\link[sp]{SpatialGridDataFrame}}-\code{\link{class}} 
#' of the \pkg{sp} package.
#' @return Returns an object of \code{\link{class}} \code{data.grid}, a list with 
#' the arguments as components.
#' @seealso \code{\link{grid.par}}, \code{\link{binning}}, \code{\link{locpol}}.
#' @export
data.grid <- function(..., grid = NULL) {
#-------------------------------------------------------------------- 
    args <- list(...)
    nargs <- length(args)
    if ( is.null(grid) ) {
        n <- dim( args[[1]] )
        if (is.null(n))
          stop("argument 'grid' (or array data) must be provided.")
          # stop("cannot derive grid parameters from data!")
        grid <- grid.par(n, min = rep(1,length(n)), max = n)
    }
    if (!inherits(grid, "grid.par"))
      stop("argument 'grid' must be of class (or extending) 'grid.par'.")
    # Let's go ...
    n <- grid$n
    index <- which(sapply(args, length) == prod(n))
    if(length(index)==0) 
      stop("no data with length equal to 'prod(grid$n)'")  
    if(length(index)!=nargs)          # "not all arguments have the same length"
      warning("some data with length not equal to 'prod(grid$n)' (ignored)") 
    # NOTA: Seguramente hai unha forma mellor de facer o seguinte...
    dimres <- if(grid$nd > 1) n else NULL  # drop dimension for 1d grid
    result <- args[index]
    seqres <- seq_along(result)
    if (is.null(names(result))) names(result) <- paste("y", seqres, sep="")	
    for (i in seqres) dim(result[[i]]) <- dimres
    # Rematar de construir o obxeto
    result$grid <- grid
    oldClass(result) <- "data.grid"
    return(result)
#--------------------------------------------------------------------
} # data.grid


