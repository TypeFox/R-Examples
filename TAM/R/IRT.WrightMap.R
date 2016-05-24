


###########################################################
# S3 method WrightMap
# IRT.WrightMap <- function(object, prob.lvl= .5 , type="PV" , ...) {
IRT.WrightMap <- function(object , ...) {
    UseMethod("IRT.WrightMap")
 }
###########################################################


############################################################
IRT.WrightMap.IRT.threshold <- function( object , label.items=NULL , ... ){
		#----------------------
		# create trait distribution
		N1 <- 20000	
		thresh1 <- object
		theta <- attr(thresh1 , "theta")
		ind <- N1 * attr(thresh1 , "prob.theta")   
		TP <- nrow(theta)
		ind <- round( ind[,1] )
		theta <- theta[ rep( 1:TP , ind ) , ]
		#--------------------------
		# input for WrightMap function
		thresh0 <- as.matrix(thresh1)
		class(thresh0) <- NULL
		attr(thresh0,"prob.theta") <- attr(thresh0,"theta") <- NULL
		if ( is.null(label.items)){
			label.items <- rownames(thresh1)
							} 
		#--- create WrightMap							
        res <- WrightMap::wrightMap( thetas= theta , thresholds= thresh0 ,
				    label.items = label.items , ...)
	    invisible(res)
				}
###############################################################		



###########################################################
# Wright map for TAM models
IRT.WrightMap.tam.mml <- function( object , prob.lvl=.5 , 
		    type="PV" , ... ){
  
    #*****************************
	# extract dimensionality
	ndim <- object$ndim

    #*****************************
    # compute thresholds
	thresh <- tam.threshold( object , prob.lvl=prob.lvl )
    
	#******************************
    # person parameter estimates

	#--- WLE
	if (type=="WLE"){
	    pers.estimates <- tam.wle( object )
		if (ndim==1){
		   pers.estimates <- pers.estimates$theta	
					}
		if (ndim>1){			
			pers.estimates <- pers.estimates[ , paste0("theta.Dim0" , 1:ndim) ]		
					}
					}	
	#--- PV
	if (type=="PV"){
        pers.estimates <- WrightMap.sim.PV( object , ndim=ndim )
					}					
    #--- Population					
	if (type=="Pop"){
	  N1 <- 10000
	   if (ndim==1){
	       pers.estimates <- rep( object$theta[,1] , round( object$pi.k * N1 ) )		
					}
	   if (ndim>1){
	    	x1 <- round( object$pi.k * N1)
			theta <- object$theta
			TP <- nrow(theta)
			pers.estimates <- theta[ rep( 1:TP , x1 ) , ]
					}		   
					}
					
	#******************************
	# draw Wright Map
	res <- WrightMap::wrightMap( pers.estimates , thresh , ... )
    invisible(res)	
		}
###########################################################


#############################################################
# Wright maps for objects of class tamaan
IRT.WrightMap.tamaan <- function( object , prob.lvl=.5 , 
		      type="PV" , ... ){
		
		method <- object$tamaanify$method
				
		plot_wm <- 1 * ( method %in% c("tam.mml" , "tam.mml.2pl" ) )
				
		if ( plot_wm == 0){
			stop("Wright Map cannot be automatically created for this model!")
						}
						
		res <- IRT.WrightMap.tam.mml( object , prob.lvl=prob.lvl , type=type , ... )		
		invisible(res)
		}
###################################################################
		
###########################################################
# simulate person parameter estimates according PVs
WrightMap.sim.PV <- function( object , ndim ){
     person <- object$person
	 N <- nrow(person)
	 if (ndim==1){
	     pers.est <- stats::rnorm( N , mean=person$EAP , sd=person$SD.EAP )
				  }
	 if (ndim>1){
        pers.est <- matrix( 0 , nrow=N , ncol=ndim)
		for (dd in 1:ndim){
			pers.est[,dd] <- stats::rnorm( N , mean= person[,paste0("EAP.Dim",dd)] , 
							sd= person[,paste0("SD.EAP.Dim",dd)] )
						}
				}
	 return(pers.est)
			}
############################################################		
