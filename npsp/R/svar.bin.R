#--------------------------------------------------------------------
#   svar.bin.R (npsp package)
#--------------------------------------------------------------------
#   svar.bin        S3 class and methods
#       svar.bin.default(x, y, maxlag, nlags, minlag, estimator, ...)   
#   svariso(x, y, maxlag, nlags, minlag, estimator, ...)   
#
#   (c) R. Fernandez-Casal         Last revision: Aug 2012
#--------------------------------------------------------------------
# PENDENTE:
#   - S3 generic ?
#--------------------------------------------------------------------


#--------------------------------------------------------------------
# svar.bin(x, ...)
#--------------------------------------------------------------------
#' Linear binning of semivariances
#' 
#' Creates a \code{svar.bin} (binned semivar. + grid parameters) object with 
#' linearly binned semivariances.
#'
#' @aliases svar.bin-class
#' @param  x 	object used to select a method. Usually a matrix with the 
#' coordinates of the data locations (columns correspond with dimensions and 
#' rows with data).
#' @param  ... 	further arguments passed to or from other methods.
#' @details  Currently, only isotropic semivariogram estimation is supported.
#' 
#' If parameter \code{nlags} is not specified is set to \code{101}.
#' @return Returns an S3 object of class \code{svar.bin} (extends \code{\link{bin.data}}), 
#'    a \code{\link{data.grid}} object with the following 4 components:
#' \item{biny}{array (dimension \code{nlags}) with the binned semivariances. }
#' \item{binw}{array (dimension \code{nlags}) with the bin counts (weights).}
#' \item{grid}{a \code{\link{grid.par}}-\code{\link{class}} object with the grid parameters.}
#' \item{data}{a list with 3 components:
#' \itemize{
#'    \item{\code{x} argument \code{x}.}
#'    \item{\code{y} argument \code{y}.}
#'    \item{\code{med} (weighted) mean of the (binned) semivariances.}
#' }}
#' \item{svar}{a list of 2 components:
#' \itemize{
#'    \item{\code{type} character, type of estimation (e.g. "isotropic").}
#'    \item{\code{estimator} character, estimator name (e.g. "classical").}
#' }}
#' @seealso \code{\link{np.svariso}}, \code{\link{np.svar}}, 
#' \code{\link{data.grid}}, \code{\link{binning}}, \code{\link{locpol}}.
#' @export
svar.bin <- function(x, ...) UseMethod("svar.bin")
# S3 generic function svar.bin


#--------------------------------------------------------------------
#' @rdname svar.bin
#' @aliases svar.bin.default iso.svar svariso
#' @inheritParams np.svar.default
#' @param estimator character, estimator name (e.g. "classical"). See "Details" below.
#' @method svar.bin default
#' @export
svar.bin.default <- function(x, y, maxlag = NULL, nlags = NULL, minlag = maxlag/nlags, 
                    estimator = c("classical", "modulus"), ...) {    
# Returns an S3 object of class "svar.bin" (extends "bin.data")
# Interface to the fortran routine "set_bin"
#
#   Devuelve la rejilla binning (lineal) para la estimación np de un semivariograma isotrópico
#   Se puede emplear para estimación clásica/robusta
#--------------------------------------------------------------------
    y <- as.numeric(y)
    ny <- length(y)                       # number of data
    x <- as.matrix(x)
    if ( !identical(ny, nrow(x)) )
      stop("arguments 'y' and 'x' do not have the same length")
    # Remove missing values 
    ok <- complete.cases(x, y) # observations having no missing values across x and y
    if (any(!ok)) {
        warning("missing values removed")
        x <- x[ok,]
        y <- y[ok]
        ny <- length(y)
    }    
    nd <- ncol(x)                         # number of dimensions
    if (is.null(maxlag)) 
        maxlag <- 0.55*sqrt(sum(diff(apply(x, 2, range))^2)) # 55% of largest lag
    if (is.null(nlags)) nlags <- 101      # dimension of the binning grid
    estimator <- match.arg(estimator)
    itipo <- ifelse(estimator == "modulus", 2, 0)
    # Let's go FORTRAN!
    # subroutine svar_iso_bin(nd, x, ny, y, nlags, minlag, maxlag, itipo,
    #                           bin_lag, bin_med, bin_y, bin_w)
    # itipo   = Tipo de estimador a calcular
    #       0 = promedio de las diferencias al cuadrado 
    #           (equivalente al estimador clásico)
    #       2 = reescalado del promedio de las diferencias absolutas  
    #           (equivalente al estimador robusto)
    ret <-.Fortran("svar_iso_bin", nd = as.integer(nd), x = as.double(t(x)), 
                  ny = as.integer(ny), y = as.double(y), 
                  nlags = as.integer(nlags), minlag = as.double(minlag), 
                  maxlag = as.double(maxlag), itipo = as.integer(itipo),
                  lag = double(1), med = double(1), biny = double(nlags), binw = double(nlags))                                    
    is.na(ret$biny) <- ret$binw == 0      # biny[binw == 0] <- NA
    result <- with( ret,
              data.grid(biny = biny, binw = binw,
              grid = grid.par(n = nlags, min = minlag, lag = lag, dimnames = "h")) )
    result$data <- list(x = x, y = y, med = ret$med)
    result$svar <- list(type = "isotropic", estimator = estimator)
    oldClass(result) <- c("svar.bin", "bin.data", "bin.den", "data.grid")
    return(result)
#--------------------------------------------------------------------
} # svariso, svar.bin.default



#--------------------------------------------------------------------
#' @rdname svar.bin
#' @export
svariso <- svar.bin.default
#--------------------------------------------------------------------


