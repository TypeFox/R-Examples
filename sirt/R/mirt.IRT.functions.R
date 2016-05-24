
###############################################
# mirt IRT functions
###############################################


########################################################
# likelihood
IRT.likelihood.mirt <- function( object , ... ){    
	object0 <- object
	object <- mirt.wrapper.posterior(object)	
    ll <- object$f.yi.qk
    attr(ll,"theta") <- object$theta.k
	attr(ll,"prob.theta") <- object$pi.k
	attr(ll,"G") <- 1
    return(ll)
        }
IRT.likelihood.ConfirmatoryClass <- IRT.likelihood.mirt
IRT.likelihood.ExploratoryClass <- IRT.likelihood.mirt
IRT.likelihood.SingleGroupClass <- IRT.likelihood.mirt
##########################################################		


########################################################
# posterior
IRT.posterior.mirt <- function( object , ... ){    
	object0 <- object
	object <- mirt.wrapper.posterior(object)	
    ll <- object$f.qk.yi
    attr(ll,"theta") <- object$theta.k
	attr(ll,"prob.theta") <- object$pi.k
	attr(ll,"G") <- 1
    return(ll)
        }
IRT.posterior.ConfirmatoryClass <- IRT.posterior.mirt
IRT.posterior.ExploratoryClass <- IRT.posterior.mirt
IRT.posterior.SingleGroupClass <- IRT.posterior.mirt
##########################################################			


########################################################
# irfprob
IRT.irfprob.mirt <- function( object , ... ){    
	object0 <- object
	object <- mirt.wrapper.posterior(object)
    ll <- object$probs
    attr(ll,"theta") <- object$theta.k
	attr(ll,"prob.theta") <- object$pi.k
	attr(ll,"G") <- 1
    return(ll)	
        }
IRT.irfprob.ConfirmatoryClass <- IRT.irfprob.mirt
IRT.irfprob.ExploratoryClass <- IRT.irfprob.mirt
IRT.irfprob.SingleGroupClass <- IRT.irfprob.mirt
##########################################################		