
###########################################################
###########################################################
# IRT.expectedCounts
###########################################################
###########################################################

###########################################################
# object of class tam.mml
IRT.expectedCounts.tam <- function( object , ... ){    
	ll <- aperm( object$n.ik , c(2,3,1,4) )
    attr(ll,"theta") <- object$theta
	attr(ll,"prob.theta") <- object$pi.k
	attr(ll,"G") <- object$G
	dimnames(ll)[[1]] <- colnames(object$resp)
    return(ll)	
    return(ll)
        }
IRT.expectedCounts.tam.mml <- IRT.expectedCounts.tam	
###########################################################

###########################################################
# object of class tam.mml
IRT.expectedCounts.tam.mml.3pl <- function( object , ... ){    
	ll <- aperm( object$n.ik , c(2,3,1,4) )
	dimnames(ll)[[1]] <- colnames(object$resp)
    attr(ll,"theta") <- object$theta
	attr(ll,"prob.theta") <- object$pi.k
	res <- list( "delta" = object$delta , 
	             "delta.designmatrix" = object$delta.designmatrix )
	attr(ll,"skillspace") <- res	
	attr(ll,"G") <- object$G
    return(ll)
        }
###########################################################

###########################################################
# objects of class tamaan
IRT.expectedCounts.tamaan <- function( object , ... ){
	if (object$tamaanify$method %in% c( "tam.mml" , "tam.mml.2pl")  ){
			res0 <- IRT.expectedCounts.tam( object , ... )			
			}
	if (object$tamaanify$method == "tam.mml.3pl"){
			res0 <- IRT.expectedCounts.tam.mml.3pl( object , ... )			
							}
	return(res0)	
			}
###################################################################	