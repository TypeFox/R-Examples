
###########################################################
# extracts expected counts
IRT.expectedCounts <- function(object, ...) {
    UseMethod("IRT.expectedCounts")
}
###########################################################

###########################################################
# object of class gdm
IRT.expectedCounts.gdm <- function( object , ... ){    
	ll <- aperm( object$n.ik , c(2,3,1,4) )
    attr(ll,"theta") <- object$theta.k
	attr(ll,"prob.theta") <- object$pi.k
	attr(ll,"G") <- object$G
    return(ll)
        }
###########################################################

###########################################################
# object of class din
IRT.expectedCounts.din <- function( object , ... ){
	Ilj <- object$I.lj
	D1 <- dim(Ilj)
	ll <- array( 0 , dim=c( D1[1] , 2 , D1[2] , 1) )
	ll[,2,,1] <- object$R.lj
	ll[,1,,1] <- object$I.lj - object$R.lj		
    attr(ll,"theta") <- object$attribute.patt.splitted
	attr(ll,"prob.theta") <- object$attribute.patt$class.prob
	attr(ll,"G") <- 1
    return(ll)
        }
###########################################################

###########################################################
# object of class gdina
IRT.expectedCounts.gdina <- function( object , ... ){    
	G <- object$G	
	Ilj <- object$control$I.lj
	D1 <- dim(Ilj)
	ll <- array( 0 , dim=c( D1[1] , 2 , D1[2] , G ) )
	if (G==1){
	   ll[,2,,1] <- object$control$R.lj
	   ll[,1,,1] <- Ilj - object$control$R.lj	
	          }
	if (G>1){
	   ll[,2,,] <- object$control$R.lj.gg
	   ll[,1,,] <- object$control$I.lj.gg - object$control$R.lj.gg	
	          }		
    attr(ll,"theta") <- object$attribute.patt.splitted
	attr(ll,"prob.theta") <- object$attribute.patt[ , 1:object$G ]
	attr(ll,"G") <- object$G	
    return(ll)
        }
############################################################	

###########################################################
# object of class slca
IRT.expectedCounts.slca <- function( object , ... ){
    ll <- aperm( object$n.ik , c(2,3,1,4) )		
	res <- list( "delta" = object$delta , 
	             "delta.designmatrix" = object$delta.designmatrix )
	attr(ll,"skillspace") <- res
	attr(ll,"prob.theta") <- object$pi.k
	attr(ll,"G") <- object$G
    return(ll)
        }
############################################################

###########################################################
# object of class mcdina
IRT.expectedCounts.mcdina <- function( object , ... ){
    ll <- object$n.ik		
    attr(ll,"theta") <- object$attribute.patt.splitted
	attr(ll,"prob.theta") <- object$attribute.patt
	attr(ll,"G") <- object$G
    return(ll)
        }
############################################################	