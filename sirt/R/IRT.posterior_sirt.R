


########################################################
# posterior rasch.copula2
IRT.posterior.rasch.copula2 <- function( object , ... ){    
    ll <- object$f.qk.yi
    attr(ll,"theta") <- vec2mat.sirt( object$theta.k )
	attr(ll,"prob.theta") <- object$pi.k
	attr(ll,"G") <- 1
    return(ll)
        }
IRT.posterior.rasch.copula3 <- IRT.posterior.rasch.copula2		
##########################################################		



##########################################################
# posterior rasch.mml2
IRT.posterior.rasch.mml <- function( object , ... ){    
    ll <- object$f.qk.yi
    attr(ll,"theta") <- vec2mat.sirt( object$theta.k )
	attr(ll,"prob.theta") <- object$pi.k
	attr(ll,"G") <- object$G
    return(ll)
        }
##########################################################		
# smirt
IRT.posterior.smirt <- IRT.posterior.rasch.mml
##########################################################		
# gom
IRT.posterior.gom <- IRT.posterior.rasch.mml
##########################################################		
# rm.facets
IRT.posterior.rm.facets <- IRT.posterior.rasch.mml
##########################################################	
# rm.sdt
IRT.posterior.rm.sdt <- IRT.posterior.rasch.mml
##########################################################	
# prob.guttman
IRT.posterior.prob.guttman <- IRT.posterior.rasch.mml
##########################################################	

##########################################################
# posterior rasch.mirtlc
IRT.posterior.rasch.mirtlc <- function( object , ... ){    
    ll <- object$estep.res$f.qk.yi
    attr(ll,"theta") <- vec2mat.sirt( object$theta.k )
	attr(ll,"prob.theta") <- object$pi.k
	attr(ll,"G") <- object$G
    return(ll)
        }
##########################################################	