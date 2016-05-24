#--------------------------------------------------------------------
#   locpol.bin.R (npsp package)
#--------------------------------------------------------------------
#   locpol.bin  S3 class and methods
#       locpol.default(x, y, h, nbin, degree, drv,  hat.bin, ncv, set.NA, ...)
#       locpol.bin.data(x, h, degree, drv, hat.bin, ncv, ...) 
#       locpol.svar.bin(x, h, degree, drv, hat.bin, ncv, ...) 
#       locpol.bin.den(x, h, degree, drv, ncv, ...) 
#   locpolhcv(x, y, nbin, objective, degree, drv, ncv, cov, ...)  
#
#   (c) R. Fernandez-Casal
#   Creation date: Aug 2012     Last revision: Aug 2013
#--------------------------------------------------------------------
# PENDENTE:
#   - is.locpol.bin
#   - update.locpol.bin
#   - comprobar dimensiones parámetros: h, ncv
#   - revisar h, opción multiplo espaciado/dimensión rejilla units.h=, h(i,i) < 0 ?)
#   - DUP = FALSE en FORTRAN  (o .C y pasar interfaces a C?)
#   - opción de ncv vector
#--------------------------------------------------------------------


#--------------------------------------------------------------------
# locpol(x, ...)
#--------------------------------------------------------------------
#' Local polynomial estimation
#' 
#' Estimates a multidimensional regression function (and its first derivatives) 
#' using local polynomial kernel smoothing of linearly binned data.
#' @aliases locpol.bin-class locpol.bin
#' @param  x 	a (data) object used to select a method.
#' @param  ... 	further arguments passed to or from other methods (e.g. to \code{\link{hcv.data}}).
#' @return Returns an S3 object of class \code{locpol.bin} (locpol + bin data + grid par.). 
#' A \code{\link{bin.data}} object with the additional (some optional) 3 components:
#' \item{est}{vector or array (dimension \code{nbin}) with the local polynomial estimates. }
#' \item{locpol}{a list with 7 components:
#' \itemize{
#'    \item{\code{degree} degree of the polinomial.}
#'    \item{\code{h} bandwidth matrix.}
#'    \item{\code{rm} residual mean.}
#'    \item{\code{rss} sum of squared residuals.}
#'    \item{\code{ncv} number of cells ignored in each direction.}
#'    \item{\code{hat} (if requested) hat matrix of the binned data.}
#'    \item{\code{nrl0} (if appropriate) number of cells with data (\code{binw > 0}) 
#'    and missing estimate (\code{est == NA}).}
#' }}
#' \item{deriv}{(if requested) matrix of first derivatives.} 
# \item{deriv}{(\eqn{length(y) \times ndim}) matrix of first derivatives.} 
#' @seealso \code{\link{binning}}, \code{\link{data.grid}}, 
#' \code{\link{np.svariso}}, \code{\link{svar.bin}},
#' \code{\link{np.den}}, \code{\link{bin.den}}, \code{\link{hcv.data}}.
#' @export
locpol <- function(x, ...) UseMethod("locpol")
# S3 generic function locpol
# Returns an S3 object of class "locpol.bin" (extends "bin.data")
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# locpol.default(x, y, h = NULL, nbin = NULL, degree = 1,
#               drv = FALSE, hat.bin = FALSE, ncv = 0, set.NA = FALSE)
#--------------------------------------------------------------------
#' @rdname locpol
#' @method locpol default
#' @param  y vector of data (response variable).
#' @param  h (full) bandwidth matrix (controls the degree of smoothing). 
#' @param  nbin vector with the number of bins on each dimension.
#' @param  degree degree of the local polynomial used. Defaults to 1 (local linear estimation).
#' @param  drv logical; if \code{TRUE}, the matrix of estimated first derivatives is returned.
#' @param  hat.bin logical; if \code{TRUE}, the hat matrix of the binned data is returned.
#' @param  ncv integer; determines the number of cells leaved out in each dimension.
#' Defaults to 0 (the full data is used) and it is not normally changed by the user
#' in this setting. See "Details" below.
#' @param  set.NA logical. If \code{TRUE}, sets binning cells without data to missing.
#' @details Standard generic function with a default method (interface to the 
#' fortran routine \code{lp_raw}), in which argument \code{x} 
#' is a vector or matrix of covariates (e.g. spatial coordinates).
#'
#' If parameter \code{nbin} is not specified is set to \code{rep(25, ncol(x))}. 
#'
#' A multiplicative triweight kernel is used to compute the weights.
#' 
#' If \code{ncv > 0}, estimates are computed by leaving out cells with indexes within 
#' the intervals \eqn{[x_i - ncv + 1, x_i + ncv - 1]}, at each dimension i, where \eqn{x} 
#' denotes the index of the estimation position. \eqn{ncv = 1} corresponds with 
#' traditional cross-validation and \eqn{ncv > 1} with modified CV 
#' (see e.g. Chu and Marron, 1991, for the one dimensional case).
#' 
#' Setting \code{set.NA = TRUE} (equivalent to \code{biny[binw == 0] <- NA}) 
#' may be useful for plotting the binned averages \code{$biny}
#' (the hat matrix should be handled with care).
#' @references
#' Chu, C.K. and Marron, J.S. (1991) Comparison of Two Bandwidth Selectors
#'   with Dependent Errors. \emph{The Annals of Statistics}, \bold{19}, 1906-1918.
#'
#' Rupert D. and Wand M.P. (1994) Multivariate locally weighted least squares regression.
#'   \emph{The Annals of Statistics}, \bold{22}, 1346-1370.
#' @examples 
#' lp <- locpol(earthquakes[, c("lon", "lat")], earthquakes$mag, h = diag(2, 2), nbin = c(41,41))
#' simage(lp, main = "Smoothed magnitude")
#' contour(lp, add = TRUE)
#' 
#' bin <- binning(earthquakes[, c("lon", "lat")], earthquakes$mag, nbin = c(41,41))
#' lp2 <- locpol(bin, h = diag(2, 2))
#' all.equal(lp, lp2)
#' 
#' ## Alternatively:
#' ## lp <- locpolhcv(earthquakes[, c("lon", "lat")], earthquakes$mag, ncv = 4)
#' 
#' den <- locpol(as.bin.den(bin), h = diag(1, 2))
#' plot(den, log = FALSE, main = 'Estimated density')
#' @export
locpol.default <- function(x, y, h = NULL, nbin = NULL, degree = 1 + as.numeric(drv), 
          drv = FALSE, hat.bin = FALSE, ncv = 0, set.NA = FALSE, ...) {    
# Returns an S3 object of class "locpol.bin" (locpol + bin data + grid parameters)
# Interface to the fortran routine "lp_raw" (lp_module.f90)
#
# PENDENTE:
#   - valor por defecto para h
#   - comprobar parámetros: ncv
#   - revisar h, opción multiplo espaciado/dimensión rejilla 
#     (defecto en fortran con valor negativo h(1,1)??)
#   - eliminar (*) de fortran, fortran 2003
#--------------------------------------------------------------------
    y <- as.numeric(y)
    ny <- length(y)                              # number of data
    x <- as.matrix(x)
    if ( !identical(ny, nrow(x)) )
        stop("arguments 'y' and 'x' have incompatible dimensions")
    drv <- as.logical(drv)
    ncv <- as.integer(ncv)
    if (ncv < 0L)
        stop("argument 'ncv' must be a positive integer")
    degree <- as.integer(degree)
    if(!(degree %in% 0L:2L))
        stop("argument 'degree' must be 0, 1 or 2")
    if (drv && degree == 0L)
        stop("'degree' must be greater than or equal to 1 if 'drv == TRUE'")
    # Remove missing values  
    ok <- complete.cases(x, y) # observations having no missing values across x and y
    if (any(!ok)) {
        warning("missing values removed")
        x <- x[ok,]
        y <- y[ok]
        ny <- length(y)
    }    
    nd <- ncol(x)                                 # number of dimensions
    if(is.null(nbin)) nbin <- rep(25L, nd) else   # dimensions of the binning grid
      if (nd != length(nbin))
        stop("arguments 'x' and 'nbin' have incompatible dimensions")
    if(is.null(h)) { 
        stop("argument 'h' (bandwith matrix) must be provided")
        # h <- diag(1, ncol = nd, nrow = nd)      # bandwith matrix PENDENTE
    } else {
        h <- as.matrix(h)
        if (!is.numeric(h) || !all(dim(h) == nd))
          stop("bandwith 'h' is not a square numeric matrix of appropriate order")   
    }      
    nt <- prod(nbin)   
    hat.bin <- as.logical(hat.bin)
    hat <- if (hat.bin) double(nt*nt) else NA_real_
    deriv <- if (drv) rep(NA_real_, nt*nd) else NA_real_ 
    # Let's go FORTRAN!
    #  subroutine lp_raw( nd, nbin, ntbin, x, ny, y,                             
    #                     bin_min, bin_max, bin_med, bin_y, bin_w,                
    #                     h, lpe, degree, ideriv, deriv, ihat, hatlp,             
    #                     ncv, rm, rss, nrl0)
    ret <-.Fortran("lp_raw", nd = as.integer(nd), nbin = as.integer(nbin),
              nt = as.integer(nt), xt = as.double(t(x)), ny = as.integer(ny), 
              y = as.double(y), min = double(nd), max = double(nd), med = double(1),
              biny = double(nt), binw = double(nt), h = as.double(h), 
              elp = as.double(rep(NA_real_, nt)), degree = as.integer(degree),  
              ideriv = as.integer(drv), deriv = deriv, ihat = as.integer(hat.bin), 
              hat = hat, ncv = as.integer(ncv), rm = double(1), rss = double(1), 
              nrl0 = integer(1), NAOK = TRUE, PACKAGE = "npsp")
    # Construir o resultado
    if (set.NA) is.na(ret$biny) <- ret$binw == 0  # biny[binw == 0] <- NA
    result <- with( ret,
              data.grid(est = elp, biny = biny, binw = binw,
              grid = grid.par(n = nbin, min = min, max = max, 
              dimnames = dimnames(x)[[2]])) )         
    result$data <- list(x = x, y = y, med = ret$med)
    result$locpol <- with( ret, 
              list( degree = degree, h = matrix(h, nrow = nd), rm = rm, rss = rss, ncv = ncv ))
    if (hat.bin) result$locpol$hat <- matrix(ret$hat, nrow = nt)
    if (ret$nrl0 > 0) {
        warning("Not enough data in some neighborhoods ('NRL < NINDRL'): ", ret$nrl0)
        result$locpol$nrl0 <- ret$nrl0
    }    
    if (drv) {
        result$deriv <- ret$deriv
        if (nd > 1) dim(result$deriv) <- c(nbin, nd)
    }
    oldClass(result) <- c("locpol.bin", "bin.data", "bin.den", "data.grid")
    return(result)
#--------------------------------------------------------------------
} # locpol.default



#--------------------------------------------------------------------
# locpol.bin.data(x, h = NULL, degree = 1 + as.numeric(drv), drv = FALSE, 
#                 hat.bin = FALSE, ncv = 0, ...) {    
#--------------------------------------------------------------------
#' @rdname locpol
#' @method locpol bin.data
#' @export
locpol.bin.data <- function(x, h = NULL, degree = 1 + as.numeric(drv), drv = FALSE, 
                                        hat.bin = FALSE, ncv = 0, ...) {    
# Returns an object of class "locpol.bin" (from a "bin.data"-class object x)
# Interface to the fortran routine "lp_bin" (lp_module.f90)
# Cuidado x puede ser un objeto locpol.bin (con componentes opcionales)
#--------------------------------------------------------------------
    if (!inherits(x, "bin.data"))
      stop("function only works for objects of class (or extending) 'bin.data'")
    nd <- x$grid$nd
    nbin <- x$grid$n
    drv <- as.logical(drv)
    ncv <- as.integer(ncv)
    if (ncv < 0L)
        stop("argument 'ncv' must be a positive integer")
    degree <- as.integer(degree)
    if(!(degree %in% 0L:2L))
        stop("argument 'degree' must be 0, 1 or 2")
    if (drv && degree == 0L)
        stop("'degree' must be greater than or equal to 1 if 'drv == TRUE'")    
    if(is.null(h)) { 
        stop("argument 'h' (bandwith matrix) must be provided")
        # h <- diag(1, ncol = nd, nrow = nd)      # bandwith matrix PENDENTE
    } else {
        h <- as.matrix(h)
        if (!is.numeric(h) || !all(dim(h) == nd))
          stop("bandwith 'h' is not a square numeric matrix of appropriate order")   
    }      
    nt <- prod(nbin)
    hat.bin <- as.logical(hat.bin)
    hat <- if (hat.bin) double(nt*nt) else NA_real_
    deriv <- if (drv) rep(NA_real_, nt*nd) else NA_real_ 
    # Let's go FORTRAN!
    #  subroutine lp_bin( nd, nbin, ntbin, bin_min, bin_max, bin_med, bin_y, bin_w,  
    #                     h, lpe, degree, ideriv, deriv, ihat, hatlp,               
    #                     ncv, rm, rss, nrl0)
    ret <-.Fortran("lp_bin", nd = as.integer(nd), nbin = as.integer(nbin),
              nt = as.integer(nt), min = as.double(x$grid$min), max = as.double(x$grid$max), 
              med = as.double(x$data$med), biny = as.double(x$biny), 
              binw = as.double(x$binw), h = as.double(h), 
              elp = as.double(rep(NA_real_, nt)), degree = as.integer(degree),  
              ideriv = as.integer(drv), deriv = deriv, ihat = as.integer(hat.bin), 
              hat = hat, ncv = as.integer(ncv), rm = double(1), rss = double(1), 
              nrl0 = integer(1), NAOK = TRUE, PACKAGE = "npsp")
    # Construir o resultado
    dim(ret$elp) <- if(nd > 1) nbin else NULL
    # if (inherits(x, "locpol.bin")) x$est <- NULL
    x$est <- NULL
    result <- c(list(est = ret$elp), x)
    result$locpol <- with( ret, 
              list( degree = degree, h = matrix(h, nrow = nd), rm = rm, rss = rss, ncv = ncv ))
    result$locpol$hat <- if (hat.bin) matrix(ret$hat, nrow = nt) else NULL 
    result$locpol$nrl0 <- NULL      
    if (ret$nrl0 > 0) {
        warning("Not enough data in some neighborhoods ('NRL < NINDRL'): ", ret$nrl0)
        result$locpol$nrl0 <- ret$nrl0 
    }    
    if (drv) {
        result$deriv <- ret$deriv
        if (nd > 1) dim(result$deriv) <- c(nbin, nd)
    } else result$deriv <- NULL
    oldClass(result) <- c("locpol.bin", "bin.data", "bin.den", "data.grid")
    return(result)
#--------------------------------------------------------------------
} # locpol.bin.data
    


#--------------------------------------------------------------------
#' @rdname locpol
#' @method locpol svar.bin
#' @return \code{locpol.svar.bin} returns an S3 object of class \code{\link{np.svar}} 
#' (locpol semivar + bin semivar + grid par.).
#' @export
locpol.svar.bin <- function(x, h = NULL, degree = 1, drv = FALSE, 
                                        hat.bin = TRUE, ncv = 0, ...){
# @seealso \code{\link{np.svariso}}, \code{\link{svar.bin}}.
#--------------------------------------------------------------------
    result <- locpol.bin.data(x, h = h, degree = degree, drv = drv, 
                  hat.bin = hat.bin, ncv = ncv, ...)
    oldClass(result) <- c("np.svar", "svar.bin", "bin.data", "bin.den", "data.grid")
    return(result)
}



#--------------------------------------------------------------------
#' @rdname locpol
#' @method locpol bin.den
#' @return \code{locpol.bin.den} returns an S3 object of class \code{\link{np.den}} 
#' (locpol den + bin den + grid par.).
#' @export
locpol.bin.den <- function(x, h = NULL, degree = 1 + as.numeric(drv), drv = FALSE, 
                              ncv = 0, ...) {    
#--------------------------------------------------------------------
    if (!inherits(x, "bin.den"))
      stop("function only works for objects of class (or extending) 'bin.den'")
    nd <- x$grid$nd
    nbin <- x$grid$n
    drv <- as.logical(drv)
    ncv <- as.integer(ncv)
    if (ncv < 0L)
        stop("argument 'ncv' must be a positive integer")
    degree <- as.integer(degree)
    if(!(degree %in% 0L:2L))
        stop("argument 'degree' must be 0, 1 or 2")
    if (drv && degree == 0L)
        stop("'degree' must be greater than or equal to 1 if 'drv == TRUE'")    
    if(is.null(h)) { 
        stop("argument 'h' (bandwith matrix) must be provided")
        # h <- diag(1, ncol = nd, nrow = nd)      # bandwith matrix PENDENTE
    } else {
        h <- as.matrix(h)
        if (!is.numeric(h) || !all(dim(h) == nd))
          stop("bandwith 'h' is not a square numeric matrix of appropriate order")   
    }      
    nt <- prod(nbin)
    # hat.bin <- as.logical(hat.bin)
    # hat <- if (hat.bin) double(nt*nt) else NA_real_
    hat <- NA_real_
    deriv <- if (drv) rep(NA_real_, nt*nd) else NA_real_ 
    # Let's go FORTRAN!
    #  subroutine lp_data_grid( nd, nbin, ntbin, bin_min, bin_max, bin_med, bin_y,
    #                       h, lpe, degree, ideriv, deriv, ihat, hatlp,
    #                       ncv, rm, rss, nrl0)
    ret <-.Fortran("lp_data_grid", nd = as.integer(nd), nbin = as.integer(nbin),
              nt = as.integer(nt), min = as.double(x$grid$min), max = as.double(x$grid$max), 
              med = as.double(0), biny = as.double(x$binw/(sum(x$binw)*prod(x$grid$lag))), 
              h = as.double(h), elp = double(nt), degree = as.integer(degree),  
              ideriv = as.integer(drv), deriv = deriv, ihat = as.integer(0L), 
              hat = hat, ncv = as.integer(ncv), rm = double(1), rss = double(1), 
              nrl0 = integer(1), NAOK = TRUE, PACKAGE = "npsp")
    # NOTA: se podría evitar el cálculo de sum(x$binw) (ojo con svar.bin)
    # Construir o resultado
    ret$elp[ret$elp < 0] <- 0
    dim(ret$elp) <- if(nd > 1) nbin else NULL   
    x$est <- NULL
    result <- c(list(est = ret$elp), x)
    result$locpol <- with( ret, 
              list( degree = degree, h = matrix(h, nrow = nd), rm = rm, rss = rss, ncv = ncv ))
    # result$locpol$hat <- if (hat.bin) matrix(ret$hat, nrow = nt) else NULL 
    result$locpol$nrl0 <- NULL      
    # if (ret$nrl0 > 0) {
    #     warning("Not enough data in some neighborhoods ('NRL < NINDRL'): ", ret$nrl0)
    #     result$locpol$nrl0 <- ret$nrl0 
    # }    
    if (drv) {
        result$deriv <- ret$deriv
        if (nd > 1) dim(result$deriv) <- c(nbin, nd)
    } else result$deriv <- NULL
    oldClass(result) <- c("np.den", "bin.den", "data.grid")
    return(result)
#--------------------------------------------------------------------
} # locpol.bin.den



#--------------------------------------------------------------------
# locpolhcv(x, y, nbin = NULL, objective = c("CV", "GCV", "MASE"),  
#           degree = 1 + as.numeric(drv), drv = FALSE,
#           ncv = ifelse(objective == "GCV", 0, 1) , cov = NULL, ...)  
#--------------------------------------------------------------------
#' @rdname locpol  
#' @inheritParams hcv.data
#' @details  \code{locpolhcv} calls \code{\link{hcv.data}} to obtain an "optimal" 
#' bandwith (additional arguments \code{...} are passed to this function). 
#' Argument \code{ncv} is only used here at the bandwith
#' selection stage (estimation is done with all the data).
#' @export
locpolhcv <- function(x, y, nbin = NULL, objective = c("CV", "GCV", "MASE"),  
                      degree = 1 + as.numeric(drv), drv = FALSE,
                      hat.bin = FALSE, set.NA = FALSE, 
                      ncv = ifelse(objective == "GCV", 0, 1), cov = NULL, ...) { 
#--------------------------------------------------------------------
    objective <- match.arg(objective)
    bin <- binning(x, y, nbin = nbin, set.NA = set.NA)
    if(is.null(cov))
        hopt <- h.cv.bin.data(bin, objective = objective, degree = degree, ncv = ncv, 
                  cov.bin = NULL, ...)$h
    else 
        hopt <- hcv.data(bin, objective = objective, degree = degree, ncv = ncv, 
                  cov = cov, ...)$h       
    return(locpol(bin, h = hopt, degree = degree, drv = drv, ncv = 0, 
        hat.bin = hat.bin))
}                  



