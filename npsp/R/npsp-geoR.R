#--------------------------------------------------------------------
#   npsp-geoR.R (npsp package)
#   trying to KISS (keep it small and simple) geoR
#--------------------------------------------------------------------
#   as.variogram()  S3 generic 
#       as.variogram.svar.bin(x, ...)
#       as.variogram.np.svar(x, ...)
#   as.variomodel()  S3 generic 
#       as.variomodel.svarmod(m, ...)
#
#   (c) R. Fernandez-Casal         Last revision: Apr Aug 2013
#--------------------------------------------------------------------
# PENDENTE:
#   - @examples
#--------------------------------------------------------------------


#' @name npsp-geoR
#' @title Interface to package "geoR"
#' @description Utilities to interact with the \pkg{geoR} package. 
NULL


#' @rdname npsp-geoR
#' @details \code{as.variogram} tries to convert a semivariogram estimate \eqn{\hat{\gamma}(h_i)} 
#' to an object of the (not fully documented) \pkg{geoR}-class \code{variogram} 
#' (see e.g. \code{\link[geoR]{variog}}).
#' @aliases variogram
#' @param  x  semivariogram estimate (e.g. \code{\link{svar.bin}} or \code{\link{np.svar}} object).
#  @param  ... further arguments passed to or from other methods.
#' @seealso \code{\link[geoR]{variog}}, \code{\link[geoR]{variofit}}, \code{\link{variomodel}}, 
#' \code{\link{svar.bin}}, \code{\link{np.svar}}.
#' @export
as.variogram <- function(x, ...) UseMethod("as.variogram")
#--------------------------------------------------------------------


#--------------------------------------------------------------------
#' @rdname npsp-geoR
#' @method as.variogram svar.bin
#' @export
as.variogram.svar.bin <- function(x, ...) {
# Convierte el resultado de \code{\link{svariso}} en un objeto
#   "compatible" con la clase "variogram" de geoR.
#
#   PENDENTE:
#     - verificar missing values
#--------------------------------------------------------------------
    if (!inherits(x, "svar.bin"))
      stop("function only works for objects of class (or extending) 'svar.bin'.")
    if (x$svar$type != "isotropic")
      stop("geoR variogram-class only accepts isotropic variograms.")
    result <- list()
    result$u <- as.numeric(coords(x))
    result$v <- x$biny
    result$n <- x$binw # nº de aportaciones
    result$trend <- "cte"
    result$beta.ols <- x$data$med
    result$max.dist <- result$u[x$grid$n]
    result$output.type <- "bin"
    oldClass(result) <- "variogram"
    return(result)
} # as.variogram.svar.bin


#--------------------------------------------------------------------
#' @rdname npsp-geoR
#' @method as.variogram np.svar
#' @export
as.variogram.np.svar <- function(x, ...) {
# Convierte el resultado de \code{\link{np.svariso}} en un objeto
#   "compatible" con la clase "variogram" de geoR.
#
#   PENDENTE:
#     - verificar missing values
#--------------------------------------------------------------------
    if (!inherits(x, "np.svar"))
      stop("function only works for objects of class (or extending) 'np.svar'.")
    result <- as.variogram.svar.bin(x, ...)
    result$v <- x$est
    if (!is.null(x$locpol$hat)) {
        # Aproximación varianza estilo Cressie (suponiendo independencia) para ajuste wls
        # PENDIENTE: ESCRIBIR/REVISAR ESTAS CUENTAS
        result$n <- with(x, rowSums(locpol$hat^2 / matrix(binw, nrow=grid$n, ncol=grid$n, byrow=TRUE)))
        result$n <- 1/result$n    # nº equivalente de aportaciones
    } # else result$n <- x$binw # nº de aportaciones
    return(result)
} # as.variogram.svar.bin



#--------------------------------------------------------------------
# as.variomodel(m, ...)
#--------------------------------------------------------------------
#' @rdname npsp-geoR
#' @aliases variomodel
#' @details \code{as.variomodel} tries to convert a semivariogram model \eqn{\gamma(pars; h)} 
#' to an object of the \pkg{geoR}-class \code{variomodel} 
#' (see e.g. \code{\link[geoR]{variofit}}).
#' @param  m  variogram model (e.g. \code{\link{svarmod}} object).
#' @param  ... further arguments passed to or from other methods.
#' @export
as.variomodel <- function(m, ...) UseMethod("as.variomodel")
#   Convierte un modelo de variograma en un objeto
#   (more or less) "compatible" con la clase "variomodel" de geoR.
#--------------------------------------------------------------------



#--------------------------------------------------------------------
#' @rdname npsp-geoR
#' @method as.variomodel svarmod
#' @export
as.variomodel.svarmod <- function(m, ...){
# Convierte el resultado de \code{\link{svarmod}} en un objeto
# "compatible" con la clase "variomodel" de geoR.
# @seealso \code{\link[geoR]{practicalRange}}
#--------------------------------------------------------------------
    if (!inherits(m, "svarmod"))
        stop("argument 'm' must be of class (or extending) 'svarmod'.")
    tmp <- match(m$model, svarmodels("isotropic"), nomatch = 0)
    if (m$type != "isotropic" || !tmp || tmp > 8)
        stop("variogram model not supported by 'geoR'.")                
    res$cov.model <- m$model
    res$cov.pars <- m$par[1:2]
    res$nugget <- m$nugget  # par[3]
    res$kappa <- m$par[4]       
    oldClass(res) <- c("variomodel")
    return(res)
}    



##--------------------------------------------------------------------
## @rdname as.variomodel
## @method as.variomodel vgm
#as.variomodel.vgm <- function(x, ...) {
##   Convierte el resultado de \code{\link{vgm}} en un objeto
##   "compatible" con la clase "variomodel" de geoR.
##--------------------------------------------------------------------
#    res <- NULL 
#    models <- c("exponential", "spherical", "circular", "gaussian", 
#        "matern", "power", "nugget", "linear", "Shapiro-Botha", "extended Shapiro-Botha" ) 
#    names(models) <- c("Exp", "Sph", "Cir", "Gau", "Mat", "Pow", "Nug", "Lin",
#        "SB", "ESB")
#    if (is.null(x) || is.null(x$model)) return(models) 
#    res$cov.model <- x$model
#    #> levels(gstat::vgm()$long)
#    # [1] "Nug (nugget)"                             
#    # [2] "Exp (exponential)"                        
#    # [3] "Sph (spherical)"                          
#    # [4] "Gau (gaussian)"                           
#    # [5] "Exclass (Exponential class)"              
#    # [6] "Mat (Matern)"                             
#    # [7] "Mat (Matern, M. Stein's parameterization)"
#    # [8] "Cir (circular)"                           
#    # [9] "Lin (linear)"                             
#    #[10] "Bes (bessel)"                             
#    #[11] "Pen (pentaspherical)"                     
#    #[12] "Per (periodic)"                           
#    #[13] "Hol (hole)"                               
#    #[14] "Log (logarithmic)"                        
#    #[15] "Pow (power)"                              
#    #[16] "Spl (spline)"                             
#    #[17] "Leg (Legendre)"                           
#    #[18] "Err (Measurement error)"                  
#    #[19] "Int (Intercept)" 
#    if (!x$model %in% models) stop("variogram model not supported")
#    res$cov.pars <- with(x, c(sill, range))
#    res$nugget <- x$nugget
#    res$kappa <- x$par$dk       
#    oldClass(res) <- c("variomodel")
#    return(res)
#}
# Rutinas gstat:
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
    