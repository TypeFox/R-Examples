vgmSCV <- function (obj, digits = 4) UseMethod("vgmSCV")
#' Spatially correlated variance (SCV)
#' 
#' Compute the proportion of the variance that is spatially correlated.
#' 
#' @param obj Variogram model fitted with available function in geostatistical
#' packages such as \pkg{gstat}, \pkg{geoR}, and \pkg{georob}.
#' 
#' @param digits Integer indicating the number of decimal places to be used.
#' 
#' @return Numeric value indicating the proportion of the variance that is 
#' spatially correlated.
#' 
#' @author Alessandro Samuel-Rosa <\email{alessandrosamuelrosa@@gmail.com}>
#' @seealso \code{\link[pedometrics]{vgmLags}}
#' @aliases vgmSCV vgmSCV.variomodel vgmSCV.variogramModel vgmSCV.georob
#' @concept variogram
#' @name vgmSCV
#' @export
# FUNCTION - geoR ##############################################################
vgmSCV.variomodel <- 
  function (obj, digits = 4) {
    res <- obj$sigmasq / (obj$sigmasq + obj$tausq)
    res <- round(unname(res), digits = digits)
    return (res)
  }
# FUNCTION - gstat #############################################################
#' @export
#' @rdname vgmSCV
vgmSCV.variogramModel <- 
  function (obj, digits = 4) {
    res <- obj$psill[2] / sum(obj$psill)
    res <- round(unname(res), digits = digits)
    return (res)
  }
# FUNCTION - georob ############################################################
#' @export
#' @rdname vgmSCV
vgmSCV.georob <- 
  function (obj, digits = 4) {
    res <- obj$param[1] / sum(obj$param[1:3])
    res <- round(unname(res), digits = digits)
    return (res)
  }
