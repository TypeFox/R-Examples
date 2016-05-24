coef.WCE <- function(object, ...){ # obtains coefficients from WCE object
if  (is.null(object$covariates) == T) {
ret <- list(covariates = NULL, WCEest = .nicer(object$est))} 
		else {
ret <- list(WCEest = .nicer(object$est), covariates = object$beta.hat.covariates)}
ret
}