#--------------------------------------------------------------------
#   np.svar.R (npsp package)
#--------------------------------------------------------------------
#   np.svar         S3 class and methods
#   np.svar()       S3 generic
#       np.svar.default(x, y, h, maxlag, nlags, minlag, degree,
#                   drv, hat.bin, ncv, ...)
#       np.svar.svar.bin(x, h, degree, drv,  hat.bin, ncv, ...) 
#   np.svariso(x, y, h, maxlag, nlags, minlag, degree,
#                   drv, hat.bin, ncv, ...)   
#   np.svariso.hcv(x, y, maxlag, nlags, minlag, degree, drv, hat.bin,
#                   objective, ncv, cov.bin, ...)
#   np.svariso.corr(lp, x, h, maxlag, nlags, minlag, degree, drv, hat.bin,
#                   tol, max.iter, plot, ylim)   
#
#   (c) R. Fernandez-Casal         Last revision: Aug 2013
#--------------------------------------------------------------------
# PENDENTE:
#   - svarisohcv o final da documentación
#   - engadir aniso, 2iso, niso
#--------------------------------------------------------------------


#--------------------------------------------------------------------
#' Local polynomial estimation of the semivariogram
#'
#' Estimates a multidimensional semivariogram (and its first derivatives) 
#' using local polynomial kernel smoothing of linearly binned semivariances.
#' @aliases np.svar-class
#' @param  x 	object used to select a method. Usually a matrix with the 
#' coordinates of the data locations (columns correspond with dimensions and 
#' rows with data).
#' @param  ... 	further arguments passed to or from other methods.
#' @details  Currently, only isotropic semivariogram estimation is supported.
#' 
#' If parameter \code{nlags} is not specified is set to \code{101}.
#'
#' The computation of the hat matrix of the binned semivariances (\code{hat.bin = TRUE})
#' allows for the computation of approximated estimation variances (e.g. in \code{\link{fitsvar.sb.iso}}).
#'
#' A multiplicative triweight kernel is used to compute the weights.
#' @return Returns an S3 object of class \code{np.svar} (locpol svar + binned svar + grid par.), 
#'     extends \code{\link{svar.bin}}, with the additional (some optional) 3 components:
#' \item{est}{vector or array with the 
#'    local polynomial semivariogram estimates. }
#' \item{locpol}{a list of 6 components:
#' \itemize{
#'    \item{\code{degree} degree of the local polinomial used.}
#'    \item{\code{h} smoothing matrix.}
#'    \item{\code{rm} mean of residual semivariances.}
#'    \item{\code{rss} sum of squared residual semivariances.}
#'    \item{\code{ncv} number of cells ignored in each direction.}
#'    \item{\code{hat} (if requested) hat matrix of the binned semivariances.}
#'    \item{\code{nrl0} (if appropriate) number of cells with \code{binw > 0} 
#'    and \code{est == NA}.}
#' }}
#' \item{deriv}{(if requested) matrix of estimated first semivariogram derivatives.} 
#' @seealso \code{\link{svar.bin}}, \code{\link{data.grid}}, \code{\link{locpol}}.
#' @references
#' Fernandez Casal R., Gonzalez Manteiga W. and  Febrero Bande M. (2003) 
#' Space-time dependency modeling using general classes of flexible stationary 
#' variogram models, \emph{J. Geophys. Res.}, \bold{108}, 8779, 
#' doi:10.1029/2002JD002909.
#'
#' Garcia-Soidan P.H., Gonzalez-Manteiga W. and Febrero-Bande M. (2003) 
#' Local linear regression estimation of the variogram, 
#' \emph{Stat. Prob. Lett.}, \bold{64}, 169-179.
#' @export
np.svar <- function(x, ...) UseMethod("np.svar")
# S3 generic function
# Non parametric pilot estimation of an isotropic semivariogram
# Returns an S3 object of class "np.svar" (extends "svar.bin")
#--------------------------------------------------------------------


#--------------------------------------------------------------------
#' @rdname np.svar
#' @aliases iso.np.svar
#' @method np.svar default
#' @param  y vector of data (response variable).
#' @inheritParams locpol.default
#' @param  maxlag maximum lag. Defaults to 55\% of largest lag. 
#' @param  nlags number of lags. Defaults to 101. 
#' @param  minlag minimun lag. 
#' @param  hat.bin logical; if \code{TRUE}, the hat matrix of the binned semivariances is returned.
#' @param cov.bin covariance matrix of the binned semivariances. 
#' Defaults to identity. 
#' @export
np.svar.default <- function(x, y, h = NULL, maxlag = NULL, nlags = NULL,
                      minlag = maxlag/nlags, degree = 1,
                      drv = FALSE, hat.bin = TRUE, ncv = 0, ...) {
#   binning cells without data are set to missing. 
#   Devuelve estimador np del semivariograma y rejilla binning
#   Interfaz para la rutina de fortran "svar_iso_np"
#--------------------------------------------------------------------
    y <- as.numeric(y)
    ny <- length(y)                       # number of data
    x <- as.matrix(x)
    if ( !identical(ny, nrow(x)) )
      stop("arguments 'y' and 'x' have incompatible dimensions")
    drv <- as.logical(drv)
    degree <- as.integer(degree)
    if(!(degree %in% 0:2))
        stop("argument 'degree' must be 0, 1 or 2")
    if (drv && degree==0)
        stop("'degree' must be greater than or equal to 1 if 'drv == TRUE'")
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
    if(is.null(h)) { 
        stop("argument 'h' (bandwith) must be provided")
        # h <- 1                          # bandwith matrix PENDENTE
     } else if (!is.numeric(h) || length(h)!= 1L)
        stop("bandwith 'h' is not a numeric value")   
    hat.bin <- as.logical(hat.bin)
    hat <- if (hat.bin) double(nlags*nlags) else NA_real_
    deriv <- if (drv) rep(NA_real_, nlags) else NA_real_
    # Let's go FORTRAN!
    # subroutine svar_iso_np( nd, x, ny, y, nlags, minlag, maxlag,
    #                           bin_lag, bin_med, bin_y, bin_w,
    #                           h, lpe, degree, ideriv, deriv, ihat, hatlp,
    #                           ndelcv, rm, rss, nrl0)
    ret <-.Fortran("svar_iso_np", nd = as.integer(nd), x = as.double(t(x)),
              ny = as.integer(ny), y = as.double(y), nlags = as.integer(nlags),
              minlag = as.double(minlag), maxlag = as.double(maxlag),
              lag = double(1), med = double(1), biny = double(nlags),
              binw = double(nlags), h = as.double(h),
              elp = as.double(rep(NA_real_, nlags)), degree = as.integer(degree),
              ideriv = as.integer(drv), deriv = deriv, ihat = as.integer(hat.bin),
              hat = hat, ncv = as.integer(ncv), rm = double(1), rss = double(1),
              nrl0 = integer(1), NAOK = TRUE, PACKAGE = "npsp")
    if (ret$nrl0 > 0)
        warning("Not enough data in some neighborhoods ('NRL < NINDRL'): ", ret$nrl0)
    # Construir o resultado
    is.na(ret$biny) <- ret$binw == 0      # biny[binw == 0] <- NA
    names(ret$min) <- "h"
    result <- with( ret,
              data.grid(est = elp, biny = biny, binw = binw,
              grid = grid.par(n = nlags, min = minlag, lag = lag, dimnames = "h")) )
    result$data <- list(x = x, y = y, med = ret$med)
    result$svar <- list(type = "isotropic", estimator = "classical")
    result$locpol <- with( ret,
              list( degree = degree, h = h, rm = rm, rss = rss, ncv = ncv ))
    if (hat.bin) result$locpol$hat <- matrix(ret$hat, nrow = nlags)
    if (ret$nrl0 > 0) {
        warning("Not enough data in some neighborhoods ('NRL < NINDRL'): ", ret$nrl0)
        result$locpol$nrl0 <- ret$nrl0
    }
    if (drv) result$deriv <- ret$deriv
    oldClass(result) <- c("np.svar", "svar.bin", "bin.data", "bin.den", "data.grid")
    return(result)
#--------------------------------------------------------------------
} # svarisonp, iso.np.svar, np.svar.default


#' @rdname np.svar
#' @method np.svar svar.bin
#' @export
np.svar.svar.bin <- locpol.svar.bin


#' @rdname np.svar
#' @export
np.svariso <- np.svar.default


#--------------------------------------------------------------------
# np.svariso.hcv(x, y, maxlag = NULL, nlags = NULL, minlag = maxlag/nlags, 
#             degree = 1, drv = FALSE, hat.bin = TRUE, objective = c("CV", "GCV", "MASE"), 
#             ncv = ifelse(objective == "GCV", 0, 1) , cov.bin = NULL, ...)  
#--------------------------------------------------------------------
#' @rdname np.svar  
#' @inheritParams h.cv.bin.data 
#' @details  \code{np.svariso.hcv} calls \code{\link{h.cv}} to obtain an "optimal" 
#' bandwith (additional arguments \code{...} are passed to this function). 
#' Argument \code{ncv} is only used here at the bandwith selection stage 
#' (estimation is done with all the data).
#' @export
np.svariso.hcv <- function(x, y, maxlag = NULL, nlags = NULL, minlag = maxlag/nlags, 
                  degree = 1, drv = FALSE, hat.bin = TRUE, objective = c("CV", "GCV", "MASE"), 
                  ncv = ifelse(objective == "GCV", 0, 1) , cov.bin = NULL, ...) { 
#--------------------------------------------------------------------
    objective <- match.arg(objective)
    bin <- svariso(x, y, maxlag = maxlag, nlags = nlags, minlag = minlag, 
            estimator = "classical") 
    hopt <- h.cv.bin.data(bin, objective = objective, degree = degree, 
            ncv = ncv , cov.bin = cov.bin, ...)$h
    return(locpol(bin, h = hopt, degree = degree, drv = drv, hat.bin = hat.bin))
}           



#--------------------------------------------------------------------
#' @rdname np.svar
#' @param  lp local polynomial estimate of the trend function (object of class 
#'  \code{\link{locpol.bin}}).
#' @param tol convergence tolerance. The algorithm stops if the average of the 
#' relative squared diferences is less than \code{tol}. Defaults to 0.04. 
#' @param max.iter maximum number of iterations. Defaults to 10.
#' @param plot logical; if \code{TRUE}, the estimates obtained at each iteration are plotted.
#' @param ylim y-limits of the plot (if \code{plot == TRUE}).
## @param col colors for lines and points if \code{plot == TRUE}. 
#' @details  
#' \code{np.svariso.corr} computes a bias-corrected nonparametric semivariogram 
#' estimate using an iterative algorithm similar to that described in  
#' Fernandez-Casal and Francisco-Fernandez (2014). This procedure tries to correct
#' the bias due to the direct use of residuals (obtained in this case from a 
#' nonparametric estimation of the trend function) in semivariogram estimation.
#' @references
#' Fernandez-Casal R. and Francisco-Fernandez M. (2014) 
#' Nonparametric bias-corrected variogram estimation under non-constant trend, 
#' \emph{Stoch. Environ. Res. Ris. Assess}, \bold{28}, 1247-1259.
#' @export
np.svariso.corr <- function(lp, x = lp$data$x, h = NULL, maxlag = NULL, nlags = NULL,
                      minlag = maxlag/nlags, degree = 1, drv = FALSE, hat.bin = TRUE, 
                      tol = 0.05, max.iter = 10, plot = FALSE, 
                      ylim = c(0,2*max(svar$biny, na.rm = TRUE))) {
#--------------------------------------------------------------------
    if (!inherits(lp, "locpol.bin"))
      stop("function only works for objects of class (or extending) 'locpol.bin'")
    if (is.null(lp$locpol$hat))
      stop("'lp' must have a '$locpol$hat' component")
    
    ny <- nrow(x)
    nd <- ncol(x)                         # number of dimensions
    if (is.null(maxlag)) 
        maxlag <- 0.55*sqrt(sum(diff(apply(x, 2, range))^2)) # 55% of largest lag
    if (is.null(nlags)) nlags <- 101      # dimension of the binning grid
    if(is.null(h)) { 
        stop("argument 'h' (bandwith) must be provided")
        # h <- 1                          # bandwith matrix PENDENTE
     } else if (!is.numeric(h) || length(h)!= 1L)
        stop("bandwith 'h' is not a numeric value")   
  
    lpdat <- predict(lp, hat.data = TRUE)
    lp.resid <- lp$data$y - lpdat$y.est
    hat.trend <- lpdat$y.hat

    svar <- np.svariso(x, lp.resid, h = h, maxlag = maxlag, nlags = nlags,  
                      degree = degree, drv = drv, hat.bin = hat.bin)
    sv.lags <- coords(svar)   

    if(plot) {
        col <- rainbow(max.iter)
        plot(sv.lags, svar$biny, xlab = "distance", ylab = "semivariance", 
            ylim = ylim, col =  col[1])
        lines(sv.lags, svar$est, col =  col[1])
    }

    svar.biased <- svar$biny
    svarold <- 0
    dists <- as.vector(dist(x))   # lower triangle of the distance matrix
    for (iter in 2:max.iter) {
    # iter <- 1; iter <- iter +1

        # cov.est <- varcov(svar, coords = x, sill = 2*max(svar$est))
        cov.est <- varcov(svar, coords = x)
        cov.bias.est <- hat.trend %*% cov.est
        cov.bias.est <- cov.bias.est %*% t(hat.trend) - cov.bias.est - t(cov.bias.est)
    
        cov.bias.diag <- diag(cov.bias.est) / 2
        svar.bias <- cov.bias.diag + rep(cov.bias.diag, each = ny) - cov.bias.est
        svar.bias <- svar.bias[lower.tri(svar.bias)]

        tmp <- binning(dists, svar.bias, nbin = 2*svar$grid$n)
        tmp <- approx(coords(tmp), tmp$biny, sv.lags)$y
        svar$biny <- svar.biased - tmp
        # if (hat.bin && !drv) svar$est <- svar$locpol$hat %*% svar$biny 
        svar <- locpol(svar, h = h,  
                      degree = degree, drv = drv, hat.bin = hat.bin)
        # CUIDADO posible división por 0              
        error <- sqrt(mean((svarold/svar$est - 1)^2, na.rm = TRUE))
        if(plot) {
            lines(sv.lags, svar$est, col =  col[iter])
            points(sv.lags, svar$biny, col =  col[iter]) 
            cat('Iteration ', iter, ': ', error, '\n')
        }    
    
        if (error < tol) break
        svarold <- svar$est
    } # for (iter in 2:sv.niter) 
    if (iter == max.iter)
        warning('The maximun number of iterations has been reached.')     
    if(plot) {
        lines(sv.lags, svar$est, lwd = 2)
        points(sv.lags, svar$biny)
        legend("topleft", paste("Iteration", seq(iter)), col = c(col[seq(iter-1)], "black"), lty = 1) 
    }
    svar$svar$estimator <- "bias-corrected (residuals based)"    
    svar$svar$iter <- iter
    svar$svar$error <- error
    return(svar)
#--------------------------------------------------------------------
} # np.svariso.corr
 
