#--------------------------------------------------------------------
#   interp.R (npsp package)
#--------------------------------------------------------------------
#   interp()    S3 generic
#       interp.grid.par(object, data, newx, ...)
#       interp.data.grid(object, data.ind, newx, ...)
#       predict.locpol.bin(object, newx, hat.data, ...)
#--------------------------------------------------------------------


#--------------------------------------------------------------------
# interp(object, ...)
#--------------------------------------------------------------------
#' Fast linear interpolation of a regular grid
#'
#' Computes a linear interpolation of multidimensional regularly gridded data.
#'
#' @param object (gridded data) object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @param data vector or array of data values.
#' @param newx vector or matrix with the (irregular) locations to interpolate. 
#'    Columns correspond with dimensions and rows with data.
#' @details \code{interp} methods are interfaces to the fortran routine \code{interp_data_grid}  
#'    (in \code{grid_module.f90}).
#' @note Linear extrapolation is performed from the end nodes of the grid.
#'
#    WARNING: May fail with missing values.    
#' @return A list with two components:
#' \item{x}{interpolation locations.}
#' \item{y}{interpolated values.}
#' @seealso  \code{\link[fields]{interp.surface}}.
#' @export
interp <- function(object, ...) UseMethod("interp")
# S3 generic function interp
#--------------------------------------------------------------------


##--------------------------------------------------------------------
## @rdname interp  
## @method interp default
## @export
#interp.default <- function(object, ...)
#    stop(paste("method 'interp' is not defined for this class:", data.class(object)))
##--------------------------------------------------------------------


#--------------------------------------------------------------------
#' @rdname interp  
#' @method interp grid.par
#' @export
interp.grid.par <- function(object, data, newx, ...) {
# Linear interpolation of a multidimensional regular grid
# \code{interp} methods are interfaces to the fortran routine "interp_data_grid" (grid_module.f90)
#--------------------------------------------------------------------
#   PENDENTE: FALLA CON DATOS MISSING
#   Determinar valores newx que están dentro de cuadriculas con datos
#   Asignar NA si se pretende extrapolar      
#--------------------------------------------------------------------
    if (!inherits(object, "grid.par"))
        stop("function only works for objects of class (or extending) 'grid.par'")
    nd <- object$nd  
    n <- object$n
    nt <- prod(n)
    if(any(is.na(data)))      #  anyNA() R >= 3.1.0
        stop("'data' has missing values")     
    if(length(data) != nt) 
        stop("'length(data)' does not match the dimension of the grid 'object'")     
    newx <- as.matrix(newx)
    ny <- nrow(newx)
    if (nd != ncol(newx))
        stop("arguments 'object' and 'newx' have incompatible dimensions")    
    # Let's go FORTRAN!
    #   subroutine interp_data_grid( nd, nbin, bin_min, bin_max, ngrid, gy, newx, ny, y)
    ret <-.Fortran("interp_data_grid", nd = as.integer(object$nd), n = as.integer(n),
                  min = as.double(object$min), max = as.double(object$max),
                  nt = as.integer(nt), data = as.double(data),
                  newx = as.double(t(newx)), ny = as.integer(ny), y = double(ny), PACKAGE = "npsp") # NAOK = FALSE)
    # Construir o resultado
    return(list(x = newx, y = ret$y))
#--------------------------------------------------------------------
} # interp.grid.par


#--------------------------------------------------------------------
#' @rdname interp  
#' @method interp data.grid
#' @param data.ind integer or character with the index or name of the data component.
#' @export
interp.data.grid <- function(object, data.ind = 1, newx, ...) {
#--------------------------------------------------------------------
    if (!inherits(object, "data.grid"))
        stop("function only works for objects of class (or extending) 'data.grid'")
    return(interp.grid.par(object$grid, object[[data.ind]], newx, ...))
#--------------------------------------------------------------------
} # interp.grid.par
    


#--------------------------------------------------------------------
#' @rdname interp  
#' @method predict locpol.bin 
#' @param hat.data logical; if \code{TRUE} (and possible), the hat matrix corresponding 
#' to the (original) data is returned. 
#' @return If \code{newx == NULL}, \code{predict.locpol.bin} returns the estimates 
#' (and optionally the hat matrix) corresponding to the data
#' (otherwise \code{interp.data.grid} is called). 
#' @details \code{predict.locpol.bin} is an interface to the fortran routine 
#' \code{predict_lp} (in \code{lp_module.f90}).
#' @note WARNING: May fail with missing values (especially if \code{object$locpol$ncv > 0}).    
#' @export
predict.locpol.bin <- function(object, newx = NULL, hat.data = FALSE, ...) {
#--------------------------------------------------------------------
    if (!inherits(object, "locpol.bin"))
      stop("function only works for objects of class (or extending) 'locpol.bin'")
#    if(any(is.na(object$est))) 
#        stop("binning estimates with missing values")
    hat.data <- as.logical(hat.data)
    if (!is.null(newx)) {
        if (hat.data) warning("argument 'hat.data' ignored ('newx != NULL')")
        return(interp.data.grid(object, data.ind = 'est', newx = newx, ...))
    }    
    nd <- object$grid$nd
    nbin <- object$grid$n
    nt <- prod(nbin)
    ny <- length(object$data$y)
    if (hat.data && is.null(object$locpol$hat))
      stop("'object' must have a '$locpol$hat' component when 'hat.data == TRUE'")
    yhat <- if (hat.data) double(ny*ny) else NA_real_
    # Let's go FORTRAN!
    # subroutine predict_locpol( nd, nbin, bin_min, bin_max, bin_med, bin_y, bin_w,
    #                       ngrid, lpe, ihat, hatlp, x, ny, lpy, hatlpy)
    ret <-.Fortran("predict_locpol", nd = as.integer(nd), nbin = as.integer(nbin),
              min = as.double(object$grid$min), max = as.double(object$grid$max),
              med = as.double(object$data$med), biny = as.double(object$biny), 
              binw = as.double(object$binw), nt = as.integer(nt), est = as.double(object$est), 
              ihat = as.integer(hat.data), hat = if (hat.data) as.double(object$locpol$hat) else NA_real_, 
              x = as.double(t(object$data$x)), ny = as.integer(ny), yest = double(ny), yhat = yhat,
              NAOK = TRUE, PACKAGE = "npsp")
    # Construir o resultado    
    return(if(hat.data) list(y.est = ret$yest, y.hat = matrix(ret$yhat, nrow = ny)) else ret$yest)
#--------------------------------------------------------------------
} # predict.locpol.bin




#--------------------------------------------------------------------
#' @rdname npsp-internals
#' @method residuals locpol.bin
#' @keywords internal
#' @export
residuals.locpol.bin <- function(object, ...) {
#--------------------------------------------------------------------
    if (!inherits(object, "locpol.bin"))
      stop("function only works for objects of class (or extending) 'locpol.bin'")
    return(object$data$y - predict(object))
}#--------------------------------------------------------------------
