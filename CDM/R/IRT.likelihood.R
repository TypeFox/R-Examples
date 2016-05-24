
###########################################################
# extracts the individual likelihood
IRT.likelihood <- function (object, ...) {
    UseMethod("IRT.likelihood")
}
###########################################################

#....................................
# How to remove attributes:
# l6 <- IRT.likelihood(mod1)
# attributes(l6) <- NULL
#....................................

###########################################################
# object of class din
IRT.likelihood.din <- function( object , ... ){
    ll <- object$like
    attr(ll,"theta") <- object$attribute.patt.splitted
	attr(ll,"prob.theta") <- object$attribute.patt$class.prob
	attr(ll,"G") <- 1
    return(ll)
        }
###########################################################

###########################################################
# object of class gdina
IRT.likelihood.gdina <- function( object , ... ){
    ll <- object$like
    attr(ll,"theta") <- object$attribute.patt.splitted
	attr(ll,"prob.theta") <- object$attribute.patt[ , 1:object$G ]
	attr(ll,"G") <- object$G	
    return(ll)
        }
############################################################		

###########################################################
# object of class mcdina
IRT.likelihood.mcdina <- function( object , ... ){
    ll <- object$like
    attr(ll,"theta") <- object$attribute.patt.splitted
	attr(ll,"prob.theta") <- object$attribute.patt
	attr(ll,"G") <- object$G
    return(ll)
        }
############################################################	

###########################################################
# object of class gdm
IRT.likelihood.gdm <- function( object , ... ){
    ll <- object$p.xi.aj
    attr(ll,"theta") <- object$theta.k
	attr(ll,"prob.theta") <- object$pi.k
	attr(ll,"G") <- object$G
    return(ll)
        }
############################################################

###########################################################
# object of class slca
IRT.likelihood.slca <- function( object , ... ){
    ll <- object$p.xi.aj
    attr(ll,"theta") <- NA
	res <- list( "delta" = object$delta , 
	             "delta.designmatrix" = object$delta.designmatrix )
	attr(ll,"skillspace") <- res
	attr(ll,"prob.theta") <- object$pi.k
	attr(ll,"G") <- object$G
    return(ll)
        }
############################################################