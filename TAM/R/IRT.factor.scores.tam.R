
###########################################################
# object of class tam
IRT.factor.scores.tam <- function( object , type="EAP" , ... ){
	if ( ! ( type %in% c("MLE","WLE","EAP") ) ){
		stop("Requested type is not supported!\n")
					}
    # EAP					
    if ( type=="EAP"){
		ll <- object$person
		ll <- ll[ , grep( "EAP" , colnames(ll) ) ]
					}
	# WLE or MLE					
	if (type %in% c("WLE","MLE")){	
	   cat("***** Estimate" , type , "*****\n")
		pers <- tam.wle( object , WLE = ( type == "WLE") , ... )
		ind1 <- grep( "theta" , colnames(pers) )
		ind2 <- grep( "error" , colnames(pers) )
		ll <- pers[ , sort(c(ind1,ind2)) ]
					}						
	attr(ll,"type") <- type
    return(ll)
        }
IRT.factor.scores.tam.mml <- IRT.factor.scores.tam		
###########################################################


###########################################################
# object of class tam.mml.3pl
IRT.factor.scores.tam.mml.3pl <- function( object , type="EAP" , ... ){
	# x1 <- c("MLE","WLE","EAP")
	# only EAPs are admissible 
	x1 <- c("EAP")
	if ( ! ( type %in% x1 ) ){
		stop("Requested type is not supported!\n")
					}
    # EAP					
    if ( type=="EAP"){
		ll <- object$person
		ll <- ll[ , grep( "EAP" , colnames(ll) ) ]
					}
#	# WLE or MLE					
#	if (type %in% c("WLE","MLE")){	
#	   cat("***** Estimate" , type , "*****\n")
#		pers <- tam.wle( object , WLE = ( type == "WLE") , ... )
#		ind1 <- grep( "theta" , colnames(pers) )
#		ind2 <- grep( "error" , colnames(pers) )
#		ll <- pers[ , sort(c(ind1,ind2)) ]
#					}						
	attr(ll,"type") <- type
    return(ll)
        }
#########################################################


###########################################################
# objects of class tamaan
IRT.factor.scores.tamaan <- function( object , type="EAP" , ... ){
	if (object$tamaanify$method %in% c( "tam.mml" , "tam.mml.2pl")  ){
			res0 <- IRT.factor.scores.tam( object , type=type , ... )			
			}
	if (object$tamaanify$method == "tam.mml.3pl"){
			res0 <- IRT.factor.scores.tam.mml.3pl( object , type=type,  ... )			
							}
	return(res0)	
			}
###################################################################			
		