#--------------------------------------------------------------------
#   bin.data.R (npsp package)
#--------------------------------------------------------------------
#   bin.data  S3 class and methods
#   binning(x, y, nbin, set.NA)
#
# PENDENTE:
#   - as.bin.data.data.grid, as.bin.data.default, ...
#   - ndim() ?
#
#   (c) R. Fernandez-Casal         Last revision: Aug 2012
#--------------------------------------------------------------------


#--------------------------------------------------------------------
# binning <- function(x, y, nbin = NULL, set.NA = FALSE) {
#--------------------------------------------------------------------
#' Linear binning
#' 
#' Discretizes the data into a regular grid (computes a binned approximation) 
#' using the multivariate linear binning technique described in Wand (1994).
#'
#' @aliases bin.data-class bin.data
# @inheritParams locpol.default
#' @param  x vector or matrix of covariates (e.g. spatial coordinates). 
#' Columns correspond with covariates (coordinate dimension) and rows with data.
#' @param  y vector of data (response variable).
#' @param  nbin vector with the number of bins on each dimension.
#' @param  set.NA logical. If \code{TRUE}, sets binning cells without data to missing.
#' @details If parameter \code{nbin} is not specified is set to \code{rep(25, ncol(x))}.
#' 
#' Setting \code{set.NA = TRUE} (equivalent to \code{biny[binw == 0] <- NA}) 
#' may be useful for plotting the binned averages \code{$biny}
#' (the hat matrix should be handled with care when using \code{\link{locpol}}).
#' @return If \code{y != NULL}, an S3 object of \code{\link{class}} \code{bin.data} 
#' (gridded binned data; extends \code{\link{bin.den}}) is returned. 
#' A \code{\link{data.grid}} object with the following 4 components:
#' \item{biny}{vector or array (dimension \code{nbin}) with the bin averages. }
#' \item{binw}{vector or array (dimension \code{nbin}) with the bin counts (weights).}
#' \item{grid}{a \code{\link{grid.par}}-\code{\link{class}} object with the grid parameters.}
#' \item{data}{a list with 3 components:
#' \itemize{
#'    \item{\code{x} argument \code{x}.}
#'    \item{\code{y} argument \code{y}.}
#'    \item{\code{med} (weighted) mean of the (binned) data.}
#' }}
#' 
#' If \code{y == NULL}, \code{\link{bin.den}} is called and a  
#' \code{\link{bin.den}}-\code{\link{class}} object is returned.
#' @seealso \code{\link{data.grid}}, \code{\link{locpol}}, \code{\link{bin.den}}, 
#' \code{\link{h.cv}}.
#' @references
#' Wand M.P. (1994) Fast Computation of Multivariate Kernel Estimators.
#'   \emph{Journal of Computational and Graphical Statistics}, \bold{3}, 433-445.
#' @examples 
#' with(earthquakes, spoints(lon, lat, mag, main = "Earthquake data"))
#'
#' bin <- binning(earthquakes[, c("lon", "lat")], earthquakes$mag, nbin = c(30,30), set.NA = TRUE)
#'
#' simage(bin, main = "Binning averages")
#' with(earthquakes, points(lon, lat, pch = 20))
#' @export
binning <- function(x, y = NULL, nbin = NULL, set.NA = FALSE) {
#--------------------------------------------------------------------
# Returns an S3 object of class "bin.data" (bin data + grid parameters)
# Interface to the fortran routine "binning"
#
# PENDENTE:
#   - Binning multivariante
#
# binning <- function(x, ...) UseMethod("binning")
# binning.default <- function(x, y, nbin = NULL, ...) {
#--------------------------------------------------------------------
    if (is.null(y)) return(bin.den(x, nbin = nbin))
    y <- as.numeric(y)
    ny <- length(y)                               # number of data
    x <- as.matrix(x)
    if ( !identical(ny, nrow(x)) )
      stop("arguments 'y' and 'x' do not have the same length.")
    # Remove missing values  
    ok <- complete.cases(x, y) # observations having no missing values across x and y
    if (any(!ok)) {
        warning("missing values removed")
        x <- x[ok,]
        y <- y[ok]
        ny <- length(y)
    }    
    nd <- ncol(x)                                 # number of dimensions
    if(is.null(nbin)) nbin <- rep(25,nd) else     # dimensions of the binning grid
      #  if (!identical(nd, length(nbin)))        # Cuidado con double e integer (nd <- 1L)
      if (nd != length(nbin))
        stop("arguments 'x' and 'nbin' have incompatible dimensions.")
    nt <- prod(nbin)
    # Let's go FORTRAN!
    #   subroutine binning( nd, nbin, x, ny, y, bin_min, bin_max, bin_med, bin_y, bin_w)
    ret <-.Fortran( "binning", nd = as.integer(nd), nbin = as.integer(nbin),
                  xt = as.double(t(x)), ny = as.integer(ny), y = as.double(y),
                  min = double(nd), max = double(nd), med = double(1),
                  biny = double(nt), binw = double(nt), PACKAGE = "npsp")
    # Construir o resultado
    if (set.NA) is.na(ret$biny) <- ret$binw == 0  # biny[binw == 0] <- NA
    result <- with( ret,
              data.grid(biny = biny, binw = binw,
              grid = grid.par(n = nbin, min = min, max = max, 
              dimnames = dimnames(x)[[2]])) )
    result$data <- list(x = x, y = y, med = ret$med)
    oldClass(result) <- c("bin.data", "bin.den", "data.grid")
    return(result)
#--------------------------------------------------------------------
} # binning.default




#--------------------------------------------------------------------
#' @rdname binning  
#' @param object (gridded data) used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @export
as.bin.data <- function(object, ...) UseMethod("as.bin.data")
# S3 generic function as.bin.den
#--------------------------------------------------------------------



#--------------------------------------------------------------------
#' @rdname binning
#' @method as.bin.data data.grid
#' @param data.ind integer or character with the index or name of the component 
#'  containing the bin averages. 
#' @param weights.ind integer or character with the index or name of the component 
#'  containing the bin counts/weights (if not specified, they are set to 
#'  \code{as.numeric( is.finite( object[[data.ind]] ))}).
#' @export
as.bin.data.data.grid <- function(object, data.ind = 1, weights.ind = NULL, ...) {
#--------------------------------------------------------------------
    if (!inherits(object, "data.grid"))
        stop("function only works for objects of class (or extending) 'data.grid'")
    y <- object[[data.ind]]
    x <- coords(object$grid)        
    binw <- if (!is.null(weights.ind))  
        object[[weights.ind]] else rep(1, prod(object$grid$n)) # pending: check missing values       
    if (!all(index <- is.finite(y))) {
        binw <- as.numeric(index) * binw
        x <- x[index, ]
        y <- y[index]
    }    
    result <- data.grid(biny = object[[data.ind]], binw = binw, grid = object$grid) 
    result$data <- list(x = x, y = y, med = mean(y))
    oldClass(result) <- c("bin.data", "bin.den", "data.grid")
    return(result)
}
