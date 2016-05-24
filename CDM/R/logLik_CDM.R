#########################################
# Log likelihood functions
#########################################
# din class
logLik.din <- function (object, ...) {
	# extract log-likelihood
	out <- object$loglike
    # number of parameters
    attr(out, "df") <- sum( object$Npars )
	# extract number of observations
    attr(out, "nobs") <- object$I
    class(out) <- "logLik"
    return(out)
}
#########################################
# gdina class
logLik.gdina <- function (object, ...) {
	# extract log-likelihood
	out <- object$loglike
    # number of parameters
    attr(out, "df") <- sum( object$Npars )
	# extract number of observations
    attr(out, "nobs") <- sum(object$N)
    class(out) <- "logLik"
    return(out)
}
###########################################
# gdm class 
logLik.gdm <- function (object, ...) {
	# extract log-likelihood
	out <- object$loglike
    # number of parameters
    attr(out, "df") <- sum( object$Npars )
	# extract number of observations
    attr(out, "nobs") <- sum(object$N)
    class(out) <- "logLik"
    return(out)
}
#############################################
#########################################
# mcdina class
logLik.mcdina <- function (object, ...) {
	# extract log-likelihood
	out <- object$loglike
    # number of parameters
    attr(out, "df") <- sum( object$Npars )
	# extract number of observations
    attr(out, "nobs") <- sum(object$I)
    class(out) <- "logLik"
    return(out)
}
#############################################
#########################################
# slca class
logLik.slca <- function (object, ...) {
	# extract log-likelihood
	out <- object$loglike
    # number of parameters
    attr(out, "df") <- sum( object$Npars )
	# extract number of observations
    attr(out, "nobs") <- sum(object$N)
    class(out) <- "logLik"
    return(out)
}
#############################################