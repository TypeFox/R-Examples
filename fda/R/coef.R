coefficients <- function(object, ...)UseMethod('coef')

coef.fd <- function(object, ...) object$coef
coefficients.fd <- function(object, ...) object$coef

coef.fdPar <- function(object, ...) object$fd$coef
coefficients.fdPar <- function(object, ...) object$fd$coef

coef.fdSmooth <- function(object, ...) object$fd$coef
coefficients.fdSmooth <- function(object, ...) object$fd$coef

coef.Taylor <- function(object, ...) object$coef
coefficients.Taylor <- function(object, ...) object$coef
