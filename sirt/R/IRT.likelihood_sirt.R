

################################################
# sirt functions with S3 methods
# ---
#
# * rasch.copula2
# * rasch.mml
# * smirt
# * rasch.mirtlc
# * gom.em
# o rm.facets
# o rm.sdt
# o noharm.sirt (to be implemented)
################################################



########################################################
# likelihood rasch.copula2
IRT.likelihood.rasch.copula2 <- function( object , ... ){    
    ll <- object$f.yi.qk
    attr(ll,"theta") <- vec2mat.sirt( object$theta.k )
	attr(ll,"prob.theta") <- object$pi.k
	attr(ll,"G") <- 1
    return(ll)
        }
IRT.likelihood.rasch.copula3 <- IRT.likelihood.rasch.copula2		
##########################################################		



##########################################################
# likelihood rasch.mml2
IRT.likelihood.rasch.mml <- function( object , ... ){    
    ll <- object$f.yi.qk
    attr(ll,"theta") <- vec2mat.sirt( object$theta.k )
	attr(ll,"prob.theta") <- object$pi.k
	attr(ll,"G") <- object$G
    return(ll)
        }
##########################################################		
# smirt
IRT.likelihood.smirt <- IRT.likelihood.rasch.mml
##########################################################		
# gom.em
IRT.likelihood.gom <- IRT.likelihood.rasch.mml
##########################################################		
# rm.facets
IRT.likelihood.rm.facets <- IRT.likelihood.rasch.mml
##########################################################	
# rm.sdt
IRT.likelihood.rm.sdt <- IRT.likelihood.rasch.mml
##########################################################
# prob.guttman
IRT.likelihood.prob.guttman <- IRT.likelihood.rasch.mml
##########################################################

##########################################################
# likelihood rasch.mirtlc
IRT.likelihood.rasch.mirtlc <- function( object , ... ){    
    ll <- object$estep.res$f.yi.qk
    attr(ll,"theta") <- vec2mat.sirt( object$theta.k )
	attr(ll,"prob.theta") <- object$pi.k
	attr(ll,"G") <- object$G
    return(ll)
        }
##########################################################		