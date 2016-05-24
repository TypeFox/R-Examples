
###########################################################
IRT.factor.scores <- function (object, ...) {
    UseMethod("IRT.factor.scores")
       }
###########################################################


###########################################################
# object of class din
IRT.factor.scores.din <- function( object , type="MLE" , ... ){
	patt1 <- object$pattern
	K <- ncol( object$q.matrix )
	N <- nrow( object$pattern )
	make.split <- FALSE
	if ( ! ( type %in% c("MLE","MAP","EAP") ) ){
		stop("Requested type is not supported!\n")
					}	
	ll <- matrix( 0 , nrow=N , ncol=K )
	if (type == "MLE"){
		class1 <- paste(patt1$mle.est)
		colnames(ll) <- paste0("MLE.skill" , 1:K)		
		make.split <- TRUE
					}
	if (type == "MAP"){
		class1 <- paste(patt1$map.est)
		make.split <- TRUE
		colnames(ll) <- paste0("MAP.skill" , 1:K)
					}
	if (type == "EAP"){
		ind <- grep( "attr" , colnames(patt1) )
		ll <- patt1[ , ind ]
		colnames(ll) <- paste0("EAP.skill" , 1:K)
					}						
	if ( make.split){
		for (kk in 1:K){	
			ll[,kk] <- as.numeric(substring( class1 , kk , kk ))
						}
					}	
	attr(ll,"type") <- type		
    return(ll)
        }
###########################################################



###########################################################
# object of class gdina
IRT.factor.scores.gdina <- function( object , type="MLE" , ... ){
	patt1 <- object$pattern
	K <- ncol( object$q.matrix )
	N <- nrow( object$pattern )
	make.split <- FALSE
	if ( ! ( type %in% c("MLE","MAP","EAP") ) ){
		stop("Requested type is not supported!\n")
					}	
	ll <- matrix( 0 , nrow=N , ncol=K )
	if (type == "MLE"){
		class1 <- paste(patt1$mle.est)
		colnames(ll) <- paste0("MLE.skill" , 1:K)
		make.split <- TRUE
					}
	if (type == "MAP"){
		class1 <- paste(patt1$map.est)
		colnames(ll) <- paste0("MAP.skill" , 1:K)
		make.split <- TRUE
					}
	if (type == "EAP"){
		ind <- grep( "attr" , colnames(patt1) )
		ll <- patt1[ , ind ]
		colnames(ll) <- paste0("EAP.skill" , 1:K)
					}						
	if ( make.split){
		for (kk in 1:K){	
			ll[,kk] <- as.numeric(substring( class1 , kk , kk ))
						}
					}	
	attr(ll,"type") <- type		
    return(ll)
        }
###########################################################


###########################################################
# object of class mcdina
IRT.factor.scores.mcdina <- function( object , type="MLE" , ... ){
	if ( ! ( type %in% c("MLE","MAP","EAP") ) ){
		stop("Requested type is not supported!\n")
					}	
	if (type == "MLE"){
		ll <- object$MLE.class
					}
	if (type == "MAP"){
		ll <- object$MLE.class
				}
	if (type == "EAP"){
		ll <- object$MLE.class
					}						
	attr(ll,"type") <- type		
    return(ll)
        }
###########################################################




###########################################################
# object of class gdm
IRT.factor.scores.gdm <- function( object , type="EAP" , ... ){
	patt1 <- object$person
	if ( ! ( type %in% c("MLE","MAP","EAP") ) ){
		stop("Requested type is not supported!\n")
					}				
	cn1 <- colnames(patt1)
	if (type == "MLE"){
		ind <- grep( "MLE" ,  cn1 )
		ll <- patt1[, ind , drop=FALSE ]
					}
	if (type == "MAP"){
		ind <- grep( "MAP" ,  cn1 )
		ll <- patt1[, ind , drop=FALSE ]
				}
	if (type == "EAP"){
	    ind <- grep( "EAP" ,  cn1 )
		F2 <- length(ind)
		ll <- patt1[ , ind , drop=FALSE ]
					}						
	attr(ll,"type") <- type		
    return(ll)
        }
###########################################################




###########################################################
# object of class slca
IRT.factor.scores.slca <- function( object , type="MLE" , ... ){
	if ( ! ( type %in% c("MLE","MAP") ) ){
		stop("Requested type is not supported!\n")
					}	
	if (type == "MLE"){
		ll <- object$MLE.class
					}
	if (type == "MAP"){
		ll <- object$MLE.class
				}
	attr(ll,"type") <- type		
    return(ll)
        }
###########################################################
