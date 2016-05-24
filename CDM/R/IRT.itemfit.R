
###########################################################
IRT.itemfit <- function (object, ...) {
    UseMethod("IRT.itemfit")
       }
###########################################################

###########################################################
# object of class din
IRT.itemfit.din <- function( object , method="RMSEA" , ... ){
	if (method=="RMSEA"){
		ll <- object$itemfit.rmsea
					}		
	attr(ll,"method") <- method		
    return(ll)
        }
###########################################################

###########################################################
# object of class gdina
IRT.itemfit.gdina <- function( object , method="RMSEA" , ... ){
	if (method=="RMSEA"){
	    ll <- object$itemfit.rmsea
					}		
	attr(ll,"method") <- method		
    return(ll)
        }
###########################################################

###########################################################
# object of class gdm
IRT.itemfit.gdm <- function( object , method="RMSEA" , ... ){
	if (method=="RMSEA"){
	    ll <- object$itemfit.rmsea$rmsea
		names(ll) <- colnames(object$data)
					}		
	attr(ll,"method") <- method		
    return(ll)
        }
###########################################################

###########################################################
# object of class slca
IRT.itemfit.slca <- function( object , method="RMSEA" , ... ){
	if (method=="RMSEA"){
		n.ik <- object$n.ik
		probs <- object$pjk
		pi.k <- object$pi.k
		ll <- itemfit.rmsea( n.ik , pi.k , probs )$rmsea	
		names(ll) <- colnames(object$data)
					}		
	attr(ll,"method") <- method		
    return(ll)
        }
###########################################################

