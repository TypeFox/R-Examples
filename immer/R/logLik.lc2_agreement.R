

###############################################################
# log-likelihood lc2_agreement
logLik.lc2_agreement <- function (object, ...) {
	# extract log-likelihood
	out <- object$loglike
    # number of parameters
    attr(out, "df") <- object$model_output$npars
	# extract number of observations
    attr(out, "nobs") <- object$nobs
    class(out) <- "logLik"
    return(out)
}
#################################################################