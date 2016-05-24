
###############################################################
# log-likelihood immer_HRM
logLik.immer_HRM <- function (object, ...) {
	# extract log-likelihood
	out <- object$like
    # number of parameters
    attr(out, "df") <- sum(object$ic$Npars)
	# extract number of observations
    attr(out, "nobs") <- object$ic$N
    class(out) <- "logLik"
    return(out)
}
#################################################################



###############################################################
# log-likelihood immer_cml
logLik.immer_cml <- function (object, ...) {
	# extract log-likelihood
	out <- object$loglike
    # number of parameters
    attr(out, "df") <- sum(object$npars)
	# extract number of observations
    attr(out, "nobs") <- object$N
    class(out) <- "logLik"
    return(out)
}
#################################################################