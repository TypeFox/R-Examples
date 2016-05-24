#--------------------------------------------------------------------
#   npsp-gstat.R (npsp package)
#   trying to KISS (keep it small and simple) gstat
#--------------------------------------------------------------------
#   as.vgm()  S3 generic
#       as.vgm.variomodel(x, ...)
#       as.vgm.svarmod(x, ...)
#       as.vgm.sb.iso(x, h, sill, ...)
#   vgm.tab.svarmod(x, h, sill, ...)
#
#   (c) R. Fernandez-Casal         Last revision: Apr Aug 2013
#--------------------------------------------------------------------
# PENDENTE:
#   - @examples
#--------------------------------------------------------------------


#' @name npsp-gstat
#' @title Interface to package "gstat"
#' @description Utilities to interact with the \pkg{gstat} package. 
NULL


#' @rdname npsp-gstat
#' @details Tries to convert a variogram object to \code{\link[gstat]{vgm}} 
#' (\code{variogramModel}-\code{\link{class}} of \pkg{gstat} package).
#' S3 generic function. 
#'
#' \code{as.vgm.variomodel} tries to convert an object of class \code{variomodel}
#' defined in \pkg{geoR} (interface to \code{\link[gstat]{as.vgm.variomodel}}
#' defined in \pkg{gstat}).
#' @param  x   variogram model object (used to select a method).
#' @param  ... further arguments passed to or from other methods.
#' @seealso
#' \code{\link[gstat]{vgm}}, \code{\link{svarmod}}.
#' @export
#--------------------------------------------------------------------
as.vgm <- function(x, ...) UseMethod("as.vgm")
# PENDENTE: @examples
#--------------------------------------------------------------------

#--------------------------------------------------------------------
#' @rdname npsp-gstat
#' @method as.vgm variomodel
#' @export
as.vgm.variomodel <- function(x, ...) {
    if (!require(gstat)) stop("'gstat' package required.")    
    return(gstat::as.vgm.variomodel(x))
}    
# PENDENTE: renombrar rutina en gstat
# gstat::as.vgm.variomodel <- function (m) {
#    model = NULL
#    if (m$cov.model == "exponential") 
#        model = "Exp"
#    else if (m$cov.model == "circular") 
#        model = "Cir"
#    else if (m$cov.model == "gaussian") 
#        model = "Gau"
#    else if (m$cov.model == "linear") 
#        stop("no correct conversion available; use power model with power 1?")
#    else if (m$cov.model == "matern") 
#        model = "Mat"
#    else if (m$cov.model == "power") 
#        model = "Pow"
#    else if (m$cov.model == "spherical") 
#        model = "Sph"
#    else if (m$cov.model == "pure.nugget") 
#        return(vgm(m$nugget + m$cov.pars[1], "Nug", 0))
#    else stop("variogram model not supported")
#    vgm(m$cov.pars[1], model, m$cov.pars[2], m$nugget, kappa = m$kappa)
#}
#
# gstat::vgm(psill, model, range, nugget, add.to, anis, kappa = 0.5, ..., covtable)
#
#--------------------------------------------------------------------


#--------------------------------------------------------------------
#' @rdname npsp-gstat
#' @method as.vgm svarmod
#' @export
as.vgm.svarmod <- function(x, ...)  as.vgm.variomodel(as.variomodel.svarmod(x))
#--------------------------------------------------------------------


#--------------------------------------------------------------------
# vgm.tab.svarmod(x, h = seq(0, x$range, length = 1000), sill = x$sill, ...)
#--------------------------------------------------------------------
#' @rdname npsp-gstat
#' @param  h vector of lags at which the covariogram is evaluated.
#' @param  sill  sill of the covariogram (or pseudo-sill).
#' @details \code{vgm.tab.svarmod} converts a \code{svarmod} object to a 
#' \code{variogramModel}-\code{\link{class}} object of type \code{"Tab"} 
#' (one-dimensional covariance table). 
#' @export
vgm.tab.svarmod <- function(x, h = seq(0, x$range, length = 1000), sill = x$sill, ...) {
#--------------------------------------------------------------------
    if (!require(gstat)) stop("'gstat' package required.")    
    if (!inherits(x, "svarmod"))
        stop("argument 'x' must be of class (or extending) 'svarmod'.")
    if (x$type != "isotropic")
        stop("'gstat' variogram model 'Tab' only accepts isotropic variograms.")        
    return(gstat::vgm(model = "Tab",  covtable = cbind(h, covar(x, h, sill = sill)) ))
}
#if (model == "Tab" && !missing(covtable)) {
#    table = as.matrix(covtable)
#    if (NCOL(table) != 2) 
#        stop("covtable should be a 2-column matrix with distance and cov.")
#    range = max(table[, 1])
#    if (min(table[, 1]) != 0) 
#        stop("the first covariance value should be at distance 0.0")
#    table = table[, 2]
#    mf = factor(c(m, "Tab"), levels = c(m, "Tab"))
#    if (!missing(add.to) || !missing(nugget)) 
#        stop("cannot add submodels or nugget to covariance Table model")
#}


#--------------------------------------------------------------------
#' @rdname npsp-gstat
#' @method as.vgm sb.iso
#' @details \code{as.vgm.sb.iso} is an alias of \code{vgm.tab.svarmod}.
#' @export
as.vgm.sb.iso <- vgm.tab.svarmod
# CUIDADO:
#     as.vgm.svarmod(x, length.table = 1000, max = x$range, sill = x$sill, ...)
#     @param  length.table  number of discretization points.
#     @param  max   maximum lag at which the covariogram is evaluated.
#--------------------------------------------------------------------
