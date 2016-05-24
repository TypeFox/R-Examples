
###########################################################
# extracts the individual posterior
IRT.posterior <- function (object, ...) {
    UseMethod("IRT.posterior")
}
###########################################################



###########################################################
# object of class din
IRT.posterior.din <- function( object , ... ){
    ll <- object$posterior
    attr(ll,"theta") <- object$attribute.patt.splitted
	attr(ll,"prob.theta") <- object$attribute.patt$class.prob
	attr(ll,"G") <- 1
    return(ll)
        }
###########################################################

###########################################################
# object of class gdina
IRT.posterior.gdina <- function( object , ... ){
    ll <- object$posterior
    attr(ll,"theta") <- object$attribute.patt.splitted
	attr(ll,"prob.theta") <- object$attribute.patt[ , 1:object$G ]
	attr(ll,"G") <- object$G	
    return(ll)
        }
############################################################		

###########################################################
# object of class mcdina
IRT.posterior.mcdina <- function( object , ... ){
    ll <- object$posterior
    attr(ll,"theta") <- object$attribute.patt.splitted
	attr(ll,"prob.theta") <- object$attribute.patt
	attr(ll,"G") <- object$G
    return(ll)
        }
############################################################	

###########################################################
# object of class gdm
IRT.posterior.gdm <- function( object , ... ){
    ll <- object$posterior
    attr(ll,"theta") <- object$theta.k
	attr(ll,"prob.theta") <- object$pi.k
	attr(ll,"G") <- object$G
    return(ll)
        }
############################################################

###########################################################
# object of class slca
IRT.posterior.slca <- function( object , ... ){
    ll <- object$posterior
    attr(ll,"theta") <- NA
	res <- list( "delta" = object$delta , 
	             "delta.designmatrix" = object$delta.designmatrix )
	attr(ll,"skillspace") <- res
	attr(ll,"prob.theta") <- object$pi.k
	attr(ll,"G") <- object$G
    return(ll)
        }
############################################################