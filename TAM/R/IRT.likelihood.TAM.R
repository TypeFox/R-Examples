

###########################################################
# likelihood
# object of class tam (and tam.mml)
IRT.likelihood.tam <- function( object , ... ){
    ll <- object$like
    attr(ll,"theta") <- object$theta
	attr(ll,"prob.theta") <- object$pi.k
	attr(ll,"G") <- object$G
	attr(ll,"pid") <- object$pid
	attr(ll,"pweights") <- object$pweights	
    return(ll)
        }
IRT.likelihood.tam.mml <- IRT.likelihood.tam 		
IRT.likelihood.tam.mml.3pl <- IRT.likelihood.tam.mml 
IRT.likelihood.tam.latreg <- IRT.likelihood.tam
###########################################################


###########################################################
# objects of class tamaan
IRT.likelihood.tamaan <- function( object , ... ){
	if (object$tamaanify$method %in% c( "tam.mml" , "tam.mml.2pl")  ){
			res0 <- IRT.likelihood.tam( object , ... )			
			}
	if (object$tamaanify$method == "tam.mml.3pl"){
			res0 <- IRT.likelihood.tam.mml.3pl( object , ... )			
							}
	return(res0)	
			}
###################################################################		

###########################################################
# posterior
# object of class tam (and tam.mml)
IRT.posterior.tam <- function( object , ... ){
    ll <- object$hwt
    attr(ll,"theta") <- object$theta
	attr(ll,"prob.theta") <- object$pi.k
	attr(ll,"G") <- object$G
	attr(ll,"pid") <- object$pid
	attr(ll,"pweights") <- object$pweights		
    return(ll)
        }
IRT.posterior.tam.mml <- IRT.posterior.tam 		
IRT.posterior.tam.mml.3pl <- IRT.posterior.tam.mml 
IRT.posterior.tam.latreg <- IRT.posterior.tam 	
###########################################################

###########################################################
# objects of class tamaan
IRT.posterior.tamaan <- function( object , ... ){
	if (object$tamaanify$method %in% c( "tam.mml" , "tam.mml.2pl")  ){
			res0 <- IRT.posterior.tam( object , ... )			
			}
	if (object$tamaanify$method == "tam.mml.3pl"){
			res0 <- IRT.posterior.tam.mml.3pl( object , ... )			
							}
	return(res0)	
			}
###################################################################	