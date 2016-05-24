

###############################################################
# log-likelihood function rasch.copula2
logLik.rasch.copula2 <- function (object, ...) {
	# extract log-likelihood
	out <- - object$ic$deviance / 2 
    # number of parameters
    attr(out, "df") <- object$ic$np
	# extract number of observations
    attr(out, "nobs") <- object$ic$n
    class(out) <- "logLik"
    return(out)
}
logLik.rasch.copula3 <- logLik.rasch.copula2 
################################################################

#####################################################
# logLik.rasch.mml
logLik.rasch.mml <- logLik.rasch.copula2
#####################################################
# smirt
logLik.smirt <- logLik.rasch.copula2
#####################################################
# rasch.mirtlc
logLik.rasch.mirtlc <- logLik.rasch.copula2
#####################################################
# gom
logLik.gom <- logLik.rasch.copula2
#####################################################
# rm.facets
logLik.rm.facets <- logLik.rasch.copula2
#####################################################
# rm.sdt
logLik.rm.sdt <- logLik.rasch.copula2
#####################################################
# prob.guttman
logLik.prob.guttman <- logLik.rasch.copula2
#####################################################