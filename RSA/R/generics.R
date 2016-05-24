#' @title Return fitted values of a RSA model
#' @description
#' Return fitted values of a RSA model
#'
#' @export
#' @method fitted RSA

#' @param object An RSA object.
#' @param ... Other parameters (currently not used)
#' @param model Model on which the fitted values are based
fitted.RSA <- function(object, ..., model="full") {
	predictRSA(object, object$data[, object$IV1], object$data[, object$IV2], model=model)
}


#' @title Return residual values of a RSA model
#' @description
#' Return residual values of a RSA model
#'
#' @export
#' @method residuals RSA
#' @aliases resid

#' @param object An RSA object.
#' @param ... Other parameters (currently not used)
#' @param model Model on which the fitted values are based
residuals.RSA <- function(object, ..., model="full") {
	object$data[, object$DV] - predictRSA(object, object$data[, object$IV1], object$data[, object$IV2], model=model)
}
resid.RSA <- residuals.RSA