#--------------------------------------------------------------------
#   h.cv.R (npsp package)
#--------------------------------------------------------------------
#   h.cv()    S3 generic
#       h.cv.bin.data(bin, objective, h.start, h.lower, h.upper, degree, 
#               ncv, cov.bin, DEalgorithm, ...)
#       h.cv.bin.den(bin, h.start, h.lower, h.upper, degree, 
#               ncv, cov.bin, DEalgorithm, ...)
#   hcv.data(bin, objective, h.start, h.lower, h.upper, degree, 
#               ncv, cov, DEalgorithm, ...)
#
# PENDENTE:
#   - Documentar diferencias entre h.cv.bin.data e hcv.data
#     hcv.data(x, y, ...)?
#   - Documentar métodos, Incluir solo referencias?
#   - h.cv e locpolhcv PODEN TER PROBLEMAS CON DATOS MISSING
#     which(is.na(lp$est)) %in% which(is.na(bin$biny))
#   - opción en binning() para obtener cov binning a partir de cov datos
#   - optimizar cálculos matriciais "GCV" e "MASE"
#
#   (c) R. Fernandez-Casal         Last revision: Jan Sep 2013
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# h.cv(bin, ...)
#--------------------------------------------------------------------
#' Cross-validation methods for bandwidth selection
#'
#' Selects the bandwidth of a local polynomial kernel (regression, density or
#' variogram) estimator using (standart or modified) CV, GCV or MASE criteria.
#'
#' @param bin object used to select a method (binned data, binned density or binned semivariogram).
#' @param ... further arguments passed to or from other methods
#'    (e.g. parameters of the optimization routine).
#' @details Currently, only diagonal windows are supported.
#' @return Returns a list containing the following 3 components:
#' \item{h}{the best (diagonal) bandwidth matrix found.} 
#' \item{value}{the value of the objective function corresponding to \code{h}.}
#' \item{objective}{the criterion used.}
#' @seealso \code{\link{locpol}}, \code{\link{locpolhcv}}, \code{\link{binning}}, 
#' \code{\link{np.svar}}.
#' @export
h.cv <- function(bin, ...) UseMethod("h.cv")
#--------------------------------------------------------------------


#--------------------------------------------------------------------
# h.cv.bin.data(bin, objective = c("CV", "GCV", "MASE"), 
#               h.start = NULL, h.lower = NULL, h.upper = NULL, degree = 1, 
#               ncv = ifelse(objective == "GCV", 0, 1), cov.bin = NULL, 
#               DEalgorithm = FALSE, ...)
#--------------------------------------------------------------------
#' @rdname h.cv  
#' @method h.cv bin.data
#' @inheritParams locpol.default
#' @param objective character; optimal criterion to be used ("CV", "GCV" or "MASE").
#' @param h.start vector; initial values for the parameters (diagonal elements) to be optimized over.  
#' If \code{DEalgorithm == FALSE} (otherwise not used), defaults to \code{(3 + ncv) * lag},
#' where \code{lag = bin$grid$lag}.  
#' @param h.lower vector; lower bounds on each parameter (diagonal elements) to be optimized.  
#' Defaults to \code{(1.5 + ncv) * bin$grid$lag}.
#' @param h.upper vector; upper bounds on each parameter (diagonal elements) to be optimized.  
#'  Defaults to \code{1.5 * dim(bin) * bin$grid$lag}.
#' @param DEalgorithm logical; if \code{TRUE}, the differential evolution optimization algorithm 
#' in package \pkg{DEoptim} is used.
#' @param  ncv integer; determines the number of cells leaved out in each dimension.
#' (0 to GCV considering all the data, \eqn{>0} to traditional or modified cross-validation).
#' See "Details" bellow.
#' @param cov.bin covariance matrix of the binned data. Defaults to identity. 
#' @param warn logical; sets the handling of warning messages  
#' (normally due to the lack of data in some neighborhoods).
#' If \code{FALSE} (the default) all warnings are ignored. 
#' @details  
#' \code{h.cv} methods use binning approximations to the objective function values. 
#' If \code{ncv > 0}, estimates are computed by leaving out binning cells with indexes within 
#' the intervals \eqn{[x_i - ncv + 1, x_i + ncv - 1]}, at each dimension i, where \eqn{x} 
#' denotes the index of the estimation position. \eqn{ncv = 1} corresponds with 
#' traditional cross-validation and \eqn{ncv > 1} with modified CV 
#' (see e.g. Chu and Marron, 1991, for the one dimensional case). 
#' For standard GCV, set \code{ncv = 0} (the full data is used).
#' For theoretical MASE, set \code{y = trend.teor}, \code{cov = cov.teor} and \code{ncv = 0}.
#' 
#' If \code{DEalgorithm == FALSE}, the \code{"L-BFGS-B"} method in \code{\link{optim}} is used.
#' 
#' @references
#' Chu, C.K. and Marron, J.S. (1991) Comparison of Two Bandwidth Selectors
#'   with Dependent Errors. \emph{The Annals of Statistics}, \bold{19}, 1906-1918.
#'
#' Francisco-Fernandez M. and Opsomer J.D. (2005) Smoothing parameter selection
#'  methods for nonparametric regression with spatially correlated errors. 
#'  \emph{Canadian Journal of Statistics}, \bold{33}, 539-558.
#' @examples 
#' bin <- binning(earthquakes[, c("lon", "lat")], earthquakes$mag, nbin = c(30,30))
#' hcv <- h.cv(bin, ncv = 2)
#' lp <- locpol(bin, h = hcv$h)
#' ## Alternatively:
#' ## lp <- locpolhcv(earthquakes[, c("lon", "lat")], earthquakes$mag, nbin = c(30,30), ncv = 2)
#' 
#' simage(lp, main = 'Smoothed magnitude')
#' contour(lp, add = TRUE)
#' with(earthquakes, points(lon, lat, pch = 20))
#' 
#' ## Density estimation
#' hden <- h.cv(as.bin.den(bin))
#' den <- np.den(bin, h = hden$h)
#' 
#' plot(den, main = 'Estimated log(density)')
#' @export
h.cv.bin.data <- function(bin, objective = c("CV", "GCV", "MASE"), 
            h.start = NULL, h.lower = NULL, h.upper = NULL, degree = 1, 
            ncv = ifelse(objective == "GCV", 0, 1), cov.bin = NULL, 
            DEalgorithm = FALSE, warn = FALSE, ...) {
#--------------------------------------------------------------------
# CRITERIO APROXIMADO A PARTIR DE BINNING
# OLLO: cov.bin = MATRIZ DE COVARIANZAS DATOS BINNING
# ... parámetros adicionales para la rutina de optimización
# Actualmente solo ventana diagonal
# h.start VECTOR con aproximación inicial
#       por defecto h.start <- (3+ncv)*lag
# h.lower y h.upper VECTORES para rangos de búsqueda
#       por defecto h.lower <- (1.5+ncv)*lag
#       por defecto h.upper <- 1.5*n*lag
# Por defecto covarianza = identidad
# Para ventana "teórica": y = trend.teor, cov.bin = cov.bin.teor, ncv = 0
# Emplear DEalgorithm para asegurar convergencia al óptimo global
#
# PENDENTE:
#   - actualmente cov.bin -> mat cov datos binning(para GCV y MASE)
#   - opción en binning() para obtener cov binning a partir de cov datos
#   - optimizar cálculos matriciais "GCV" e "MASE"
#--------------------------------------------------------------------
    if (!inherits(bin, "bin.den"))
        stop("function only works for objects of class (or extending) 'bin.den'")
    nd <- bin$grid$nd
    n <- bin$grid$n
    nt <- prod(n)
    w <- bin$binw
    sw <- sum(w)     # normalmente = length(bin$data$y) salvo en semivariogramas
    objective <- match.arg(objective)
    if(is.null(cov.bin))
        # cov.bin <- diag(nt)  
        # Select objective function
        fhopt <- switch(objective,
            CV  = function(x)
                    locpol(bin, h = diag(x, nrow = nd), degree = degree, ncv = ncv)$locpol$rss/sw ,
            GCV = function(x) {
                    lp <- locpol(bin, h = diag(x, nrow = nd), degree = degree, hat.bin = TRUE, ncv = ncv)
                    return( with(lp$locpol,
                      sw*(rss / (sw - sum(w * diag(hat)))^2) ))
                } ,
            MASE= function(x) {
                    lp <- locpol(bin, h = diag(x, nrow = nd), degree = degree, hat.bin = TRUE, ncv = ncv)
                    return( with(lp$locpol,
                    (rss + sum(w * diag(tcrossprod(hat))))/sw ))  # PENDIENTE OPTIMIZAR...
                }
        ) # switch
    else {
        p <- (d <- dim(cov.bin))[1L]
        if (!is.numeric(cov.bin) || length(d) != 2L || p != d[2L] || p != nt)   # PENDENTE
          stop("'cov.bin' must be a square matrix of order 'prod(dim(bin))'")

        if  (objective == "GCV") cov.bin <- cov2cor(cov.bin)  # correlation matrix for "GCV"
        # Select objective function
        fhopt <- switch(objective,
            CV  = function(x) {
                    lp <- locpol(bin, h = diag(x, nrow = nd), degree = degree, hat.bin = TRUE, ncv = ncv)
                    return( with(lp$locpol,
                    (rss + 2 * sum(as.numeric(w) * hat * cov.bin))/sw ))  # trace(A%*%B) = sum(A*t(B))
                } ,                    
            GCV = function(x) {
                    lp <- locpol(bin, h = diag(x, nrow=nd), degree = degree, hat.bin = TRUE, ncv = ncv)
                    return( with(lp$locpol,
                      sw * (rss / (sw - sum(as.numeric(w) * hat * cov.bin))^2) ))
                } ,
            MASE= function(x) {
                    lp <- locpol(bin, h = diag(x, nrow = nd), degree = degree, hat.bin = TRUE, ncv = ncv)
                    return( with(lp$locpol,
                    (rss + sum(w * diag(hat %*% cov.bin %*% t(hat))))/sw ))
                }
        ) # switch
    } # if(is.null(cov.bin))
    # Minimization of the objective function
    lag <- bin$grid$lag
    if(is.null(h.lower)) h.lower <- (1.5+ncv)*lag
    if(is.null(h.upper)) h.upper <- 1.5*n*lag
    if(!as.logical(warn)) {
        ops <- options(warn = -1)
        on.exit(options(ops)) 
    }
    if(!DEalgorithm) {
        if(is.null(h.start)) h.start <- (3+ncv)*lag
        res <- optim( h.start, fhopt, method = "L-BFGS-B",
                      lower = h.lower, upper = h.upper, ...)
        return(list(h = diag(res$par, nrow = nd), value = res$value, objective = objective))
    } else {
        if (!require(DEoptim)) stop("'DEalgorithm' requires 'DEoptim' package")
        res <- DEoptim::DEoptim( fhopt, lower = h.lower, upper = h.upper, ...)
        return(list(h = diag(res$optim$bestmem, nrow = nd), value = res$optim$bestval, objective = objective))
    }
#--------------------------------------------------------------------
} # h.cv.bin.data 



#--------------------------------------------------------------------
# h.cv.bin.den(bin, h.start = NULL, h.lower = NULL, h.upper = NULL, 
#            degree = 1, ncv = 1, DEalgorithm = FALSE, ...)
#--------------------------------------------------------------------
#' @rdname h.cv  
#' @method h.cv bin.den
#' @export
h.cv.bin.den <- function(bin, h.start = NULL, h.lower = NULL, h.upper = NULL, 
            degree = 1, ncv = 1, DEalgorithm = FALSE, warn = FALSE, ...) {
#--------------------------------------------------------------------
    if (!inherits(bin, "bin.den") || inherits(bin, "bin.data"))
        stop("function only works for objects of class 'bin.den'")
    return(h.cv.bin.data(bin, objective = "CV", h.start = h.start,
        h.lower = h.lower, h.upper = h.upper, degree = degree, 
        ncv = ncv, cov.bin = NULL, DEalgorithm = DEalgorithm, ...))              
}            


#--------------------------------------------------------------------
# hcv.data(bin, objective = c("CV", "GCV", "MASE"), 
#         h.start = NULL, h.lower = NULL, h.upper = NULL, degree = 1, 
#         ncv = ifelse(objective == "GCV", 0, 1), cov = NULL, 
#         DEalgorithm = FALSE, ...)
#--------------------------------------------------------------------
#' @rdname h.cv  
#' @param cov covariance matrix of the data. Defaults to identity (uncorrelated data). 
#' @details  
#' \code{hcv.data} evaluates the objective functions at the original data 
#' (combining a binning approximation to the nonparametric estimates with a linear interpolation). 
#' If \code{ncv > 1} (modified CV), a similar algorithm to that in \code{h.cv.bin.data} is used, 
#' estimates are computed by leaving out binning cells with indexes within 
#' the intervals \eqn{[x_i - ncv + 1, x_i + ncv - 1]}. 
#' @export
hcv.data <- function(bin, objective = c("CV", "GCV", "MASE"), 
            h.start = NULL, h.lower = NULL, h.upper = NULL, degree = 1, 
            ncv = ifelse(objective == "GCV", 0, 1), cov = NULL, 
            DEalgorithm = FALSE, warn = FALSE, ...) {
# CRITERIO BINNING "EXACTO"
# COV = MATRIZ DE COVARIANZAS DOS DATOS ORIXINAIS
# COIDADO CO TEMPO DE CPU SE Nº DE DATOS GRANDE
#--------------------------------------------------------------------
    if (!inherits(bin, "bin.data"))
        stop("function only works for objects of class (or extending) 'bin.data'")
    if (inherits(bin, "svar.bin"))
        stop("not supported for objects of class (or extending) 'svar.bin'")
    nd <- bin$grid$nd
    n <- bin$grid$n
    ny <- length(bin$data$y)
    w <- bin$binw
    sw <- sum(w)
    objective <- match.arg(objective)
    if(is.null(cov))
        # cov <- diag(ny)
        f.h.cv <- switch(objective,
            CV  = if (ncv==1) function(x) {
                    lp <- locpol(bin, h = diag(x, nrow = nd), degree = degree, hat.bin = TRUE, ncv = 0)
                    if(!is.null(lp$locpol$nrl0))
                        lp$est[is.na(lp$est)] <- lp$data$med # evitar problemas si ncv>0
                    lpdat <- predict(lp, hat.data = TRUE)
                    return(mean(((lpdat$y.est - lp$data$y) / (1 - diag(lpdat$y.hat)))^2))
                } else function(x) {
                    lp <- locpol(bin, h = diag(x, nrow = nd), degree = degree, ncv = ncv)
                    if(!is.null(lp$locpol$nrl0))
                        lp$est[is.na(lp$est)] <- lp$data$med # evitar problemas si ncv>0
                    y.est <- predict(lp)  # hat.data = FALSE
                    return(mean((y.est - lp$data$y)^2))
                  },
            GCV = function(x) {
                    lp <- locpol(bin, h = diag(x, nrow=nd), degree = degree, hat.bin = TRUE, ncv = ncv)
                    if(!is.null(lp$locpol$nrl0))
                        lp$est[is.na(lp$est)] <- lp$data$med # evitar problemas si ncv>0
                    lpdat <- predict(lp, hat.data = TRUE)
                    return( mean((lpdat$y.est - lp$data$y)^2) / 
                        (1 - mean(diag(lpdat$y.hat)))^2 ) 
                  },
            MASE= function(x) {
                    lp <- locpol(bin, h = diag(x, nrow = nd), degree = degree, hat.bin = TRUE, ncv = ncv)
                    if(!is.null(lp$locpol$nrl0))
                        lp$est[is.na(lp$est)] <- lp$data$med # evitar problemas si ncv>0
                    lpdat <- predict(lp, hat.data = TRUE)
                    return( (sum((lpdat$y.est - lp$data$y)^2) +  sum(lpdat$y.hat^2))/ny )   
                  }
        ) # switch        
    else {
        p <- (d <- dim(cov))[1L]
        if (!is.numeric(cov) || length(d) != 2L || p != d[2L] || p != ny)
          stop("'cov' must be a square matrix of order 'length(y)'")
        if  (objective == "GCV") cov <- cov2cor(cov)  # correlation matrix for "GCV"
        # Select objective function
        f.h.cv <- switch(objective,
            CV  = if (ncv==1) function(x) {
                    lp <- locpol(bin, h = diag(x, nrow = nd), degree = degree, hat.bin = TRUE, ncv = 0)
                    if(!is.null(lp$locpol$nrl0))
                        lp$est[is.na(lp$est)] <- lp$data$med # evitar problemas si ncv>0
                    lpdat <- predict(lp, hat.data = TRUE)
                    return( (sum( ((lpdat$y.est - lp$data$y) / (1 - diag(lpdat$y.hat)))^2) 
                            + 2 * sum(lpdat$y.hat * cov ))/ny ) # OJO: Matriz de suavizado con todos los datos        
                } else function(x) {
                    lp <- locpol(bin, h = diag(x, nrow = nd), degree = degree, hat.bin = TRUE, ncv = ncv)
                    if(!is.null(lp$locpol$nrl0))
                        lp$est[is.na(lp$est)] <- lp$data$med # evitar problemas si ncv>0
                    lpdat <- predict(lp, hat.data = TRUE)
                    return( (sum((lpdat$y.est - lp$data$y)^2) + 2 * sum(lpdat$y.hat * cov ))/ny )
                  },
            GCV = function(x) {
                    lp <- locpol(bin, h = diag(x, nrow=nd), degree = degree, hat.bin = TRUE, ncv = ncv)
                    if(!is.null(lp$locpol$nrl0))
                        lp$est[is.na(lp$est)] <- lp$data$med # evitar problemas si ncv>0
                    lpdat <- predict(lp, hat.data = TRUE)
                    return( ny * ( sum((lpdat$y.est - lp$data$y)^2) / 
                            (ny - sum(lpdat$y.hat * cov))^2) )
                  },
            MASE= function(x) {
                    lp <- locpol(bin, h = diag(x, nrow = nd), degree = degree, hat.bin = TRUE, ncv = ncv)
                    if(!is.null(lp$locpol$nrl0))
                        lp$est[is.na(lp$est)] <- lp$data$med # evitar problemas si ncv>0
                    lpdat <- predict(lp, hat.data = TRUE)
                    return( mean((lpdat$y.est - lp$data$y)^2 + 
                        diag(lpdat$y.hat %*% cov %*% t(lpdat$y.hat))) )
                  }
        ) # switch
    } # if(is.null(cov))

    # Minimization of the objective function
    lag <- bin$grid$lag
    if(is.null(h.lower)) h.lower <- (1.5+ncv)*lag
    if(is.null(h.upper)) h.upper <- 1.5*n*lag
    if(!as.logical(warn)) {
        ops <- options(warn = -1)
        on.exit(options(ops)) 
    }
    if(!DEalgorithm) {
        if(is.null(h.start)) h.start <- (3+ncv)*lag
        res <- optim( h.start, f.h.cv, method = "L-BFGS-B",
                      lower = h.lower, upper = h.upper, ...)
        return(list(h = diag(res$par, nrow = nd), value = res$value, objective = objective))
    } else {
        if (!require(DEoptim)) stop("'DEalgorithm' requires 'DEoptim' package")
        res <- DEoptim::DEoptim( f.h.cv, lower = h.lower, upper = h.upper, ...)
        return(list(h = diag(res$optim$bestmem, nrow = nd), value = res$optim$bestval, objective = objective))
    }
#--------------------------------------------------------------------
} # hcv.data







