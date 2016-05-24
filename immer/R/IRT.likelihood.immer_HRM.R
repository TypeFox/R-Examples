


###########################################################
# likelihood
IRT.likelihood.immer_HRM <- function( object , ... ){
    ll <- object$f.yi.qk
    attr(ll,"theta") <- object$theta_like
	attr(ll,"prob.theta") <- object$pi.k
	attr(ll,"G") <- 1
    return(ll)
        }
#############################################################		


###########################################################
# posterior
IRT.posterior.immer_HRM <- function( object , ... ){
    ll <- object$f.qk.yi
    attr(ll,"theta") <- object$theta_like
	attr(ll,"prob.theta") <- object$pi.k
	attr(ll,"G") <- 1
    return(ll)
        }
#############################################################	