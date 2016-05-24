
###########################################################
# extracts the individual irfprob
IRT.irfprob <- function(object, ...) {
    UseMethod("IRT.irfprob")
}
###########################################################



###########################################################
# object of class din
IRT.irfprob.din <- function( object , ... ){
    ll <- object$pjk
	dimnames(ll)[[1]] <- colnames(object$data)
    attr(ll,"theta") <- object$attribute.patt.splitted
	attr(ll,"prob.theta") <- object$attribute.patt$class.prob
	attr(ll,"G") <- 1
    return(ll)
        }
###########################################################

###########################################################
# object of class gdina
IRT.irfprob.gdina <- function( object , ... ){
    ll <- object$pjk
	dimnames(ll)[[1]] <- colnames(object$data)
    attr(ll,"theta") <- object$attribute.patt.splitted
	attr(ll,"prob.theta") <- object$attribute.patt[ , 1:object$G ]
	attr(ll,"G") <- object$G	
    return(ll)
        }
############################################################		

###########################################################
# object of class mcdina
IRT.irfprob.mcdina <- function( object , ... ){
    ll <- object$pik
	dimnames(ll)[[1]] <- colnames(object$data)
    attr(ll,"theta") <- object$attribute.patt.splitted
	attr(ll,"prob.theta") <- object$attribute.patt
	attr(ll,"G") <- object$G
    return(ll)
        }
############################################################	

###########################################################
# object of class gdm
IRT.irfprob.gdm <- function( object , ... ){
    ll <- aperm( object$pjk , c(2 , 3 , 1 ) )
	dimnames(ll)[[1]] <- colnames(object$data)
    attr(ll,"theta") <- object$theta.k
	attr(ll,"prob.theta") <- object$pi.k
	attr(ll,"G") <- object$G
    return(ll)
        }
############################################################

###########################################################
# object of class slca
IRT.irfprob.slca <- function( object , ... ){
    ll <- aperm( object$pjk , c(2 , 3 , 1 ) )	
	dimnames(ll)[[1]] <- colnames(object$data)
    attr(ll,"theta") <- NA
	res <- list( "delta" = object$delta , 
	             "delta.designmatrix" = object$delta.designmatrix )
	attr(ll,"skillspace") <- res
	attr(ll,"prob.theta") <- object$pi.k
	attr(ll,"G") <- object$G
    return(ll)
        }
############################################################