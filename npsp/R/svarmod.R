#--------------------------------------------------------------------
#   svarmod.R (npsp package)
#--------------------------------------------------------------------
#   svarmod  S3 class and methods
#       svarmod(model, type, par, nugget, sill, range)
#       svarmod.sb.iso(dk, x, z, nu, range, sill)
#   svarmodels(type)
#   svm  S3 generic
#       sv.default(x, h, ...)
#       sv.svarmod(x, h, ...)
#       sv.sb.iso(x, h, ...)
#
#   (c) R. Fernandez-Casal         Last revision: Apr Sep 2013
#--------------------------------------------------------------------
# NOTA: para modelos isotrópicos paramétricos se toma (provisionalmente) 
#       como referencia geoR
#
# PENDENTE:
#   - @examples
#   - S3 method
#   - svar.vgm
#--------------------------------------------------------------------

 
#--------------------------------------------------------------------
#   svarmod(model, type = "isotropic", par = NA, nugget = 0, sill = NA, range = NA)
#--------------------------------------------------------------------
#   @aliases mod.svar mod.svar-class svar.mod-class
#   @aliases svarmod-class
#' Define a (semi)variogram model
#' 
#' Defines a variogram model specifying the parameter values. 
#' Constructor function of the \code{svarmod}-\code{\link{class}}.
#' 
#' @param  model string indicating the variogram family (see Details below).
#' @param  type  string indicating the type of variogram, e.g. "isotropic".
#' @param  par vector of variogram parameters.
#' @param  nugget  nugget value \eqn{c_0}.
#' @param  sill  variance \eqn{\sigma^2} or sill of the variogram (NA for unbounded variograms).
#' @param  range range (practical range or scale parameter) of the variogram 
#' (NA for unbounded variograms; maybe a vector for anisotropic variograms).
#' @return
#' \code{svarmod} returns an \code{svarmod}-\code{\link{class}} object, a list 
#' with function arguments as components.
#' @note \code{svarmod} does not check the consistency of the parameter values.
#' @seealso
#' \code{\link{sv}}, \code{\link{covar}}.
#' @export
#--------------------------------------------------------------------
svarmod <- function(model, type = "isotropic", par = NA,
                      nugget = NULL, sill = NULL, range = NULL) {
# Define a (semi)variogram model
# En los modelos paramétricos:
#   names(par) <- c('psill', 'phi', 'nugget', 'kappa')
#   phi = scale parameter
# PENDIENTE: ASIGNAR nugget, sill Y range AUTOMÁTICAMENTE
#--------------------------------------------------------------------
    model <- match.arg(model, svarmodels(type))
    if (model == "pure.nugget") sill <- nugget <- par[1]
    if (is.null(nugget)) nugget <- par[3]
    if (is.null(sill)) 
        sill <- ifelse(model %in% c("power", "linear"), NA, par[1] + par[3])
    if (is.null(range)) range <- switch(model,
        spherical = par[2],   
        circular = par[2],
        exponential = 3*par[2],
        gaussian = 1.73*par[2], 
        NA # default        
    )
    result <- list(model = model, type = type, par = par,
                      nugget = nugget, sill = sill, range = range)
    oldClass(result) <- c("isotropic", "svarmod")
    return(result)
} # svarmod
# svar.mod <- svarmod


#--------------------------------------------------------------------
# svarmod.sb.iso( dk, x, z, nu, range, sill = nu) 
# Define a Shapiro-Botha (semi)variogram model
# Returns an S3 object of class \code{sb.iso} (extends \code{svarmod})
#--------------------------------------------------------------------
#' @rdname svarmod  
#' @aliases sb.iso-class
#' @param  dk dimension of the kappa function.
#' @param  x  discretization nodes.
#' @param  z  jumps (of the spectral distibution) at the discretization nodes.
#' @param  nu  parameter \eqn{\nu_0} (can be thought of as the sill).
#' @return
#' \code{svarmod.sb.iso} returns an S3 object of \code{\link{class}} \code{sb.iso} 
#' (extends \code{svarmod}) corresponding to a `nonparametric' isotropic Shapiro-Botha model.
#' @references
#' Shapiro, A. and Botha, J.D. (1991) Variogram fitting with a general class of 
#'   conditionally non-negative definite functions. \emph{Computational Statistics 
#'   and Data Analysis}, \bold{11}, 87-96. 
#' @export
svarmod.sb.iso <- function( dk, x, z, nu, range, sill = nu) {
#--------------------------------------------------------------------
    result <- svarmod(model = svarmodels()["SB"], type = "isotropic",
          par = list(dk = dk, x = x, z = z, nu = nu),
          nugget = nu - sum(z), sill = sill, range = range)
    oldClass(result) <- c("sb.iso", "isotropic", "svarmod")
    return(result)
}

#--------------------------------------------------------------------
# sb.iso <- svarmod.sb.iso
# svarmod.iso.default
# svarmodsb.iso
# svar.mod
#--------------------------------------------------------------------


#--------------------------------------------------------------------
#' @rdname svarmod  
#' @return
#' \code{svarmodels} returns a named character vector with the available models 
#' of the corresponding \code{type} 
#' (when appropriate, component values could be used as \code{cov.model} argument in \pkg{geoR} routines
#'  and component names as \code{model} argument in \pkg{gstat} routines).
#' @export
svarmodels <- function(type = "isotropic") {
# match(model, svarmodels(type), nomatch = 0)
#--------------------------------------------------------------------
    models <- switch(type,
        isotropic = c(Exp = "exponential", Sph = "spherical", Cir = "circular", 
                    Gau = "gaussian", Mat = "matern", Pow = "power", Nug = "pure.nugget", 
                    Lin = "linear", SB = "Shapiro-Botha"),
        aniso2comp = c(ESB = "extended Shapiro-Botha"),
        stop("variogram type not defined.") # default        
    )
    return(models)
}


#--------------------------------------------------------------------
#   sv(x, h, ...)
#--------------------------------------------------------------------
#' Evaluate a semivariogram model 
#' 
#' Evaluates an \code{svarmod} object \code{x} at lags \code{h} (S3 generic function).
#' 
#' @param  x  variogram model (\code{\link{svarmod}} object).
#' @param  h  vector (isotropic case) or matrix of lag values.
#' @param  ... further arguments passed to or from other methods.
#' @return
#' A vector of semivariance values \eqn{\gamma(h_i)}.
#' @seealso
#' \code{\link{covar}}
#' @export
#--------------------------------------------------------------------
sv <- function(x, h, ...) UseMethod("sv")


#--------------------------------------------------------------------
#' @rdname sv
#' @method sv default
#' @export
sv.default <- function(x, h, ...) stop("Invalid variogram object")
#--------------------------------------------------------------------

#--------------------------------------------------------------------
#' @rdname sv
#' @method sv svarmod
#' @export
sv.svarmod <- function(x, h, ...) as.vgm.svarmod(x, h)$covtable
# stop("Not defined (yet) for general variogram models")
#--------------------------------------------------------------------


#--------------------------------------------------------------------
#' @rdname sv  
#' @method sv sb.iso
#' @export
sv.sb.iso <- function(x, h, ...) {
# CUIDADO SI DIMENSIONES DE h GRANDE outer(h, x)
#--------------------------------------------------------------------
    result <- with(x$par,
        drop(nu - kappasb(outer(h, x), dk) %*% z) )
    result[h < 10 * .Machine$double.eps] <- 0
    return(result)
}



