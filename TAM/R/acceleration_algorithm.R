
##############################################################################
# acceleration
accelerate_parameters <- function( xsi_acceleration , xsi , iter , itermin=2 , ind=NULL){
	
	eps <- 1E-70
	parm_history <- xsi_acceleration$parm_history	
	if ( is.null(ind) ){
		ind <- seq( 1 , nrow(parm_history) )
						}
# cat( "Iteration .... " , iter , " ..... \n") 

	
	#*********************************************************
	# overrelaxation method according to Yu (2012)
	if ( xsi_acceleration$acceleration == "Yu"){
		  w_accel <- xsi_acceleration$w		 
   	     if ( w_accel < 0 ){ w_accel <- .01 }
				# xsi <- xsi + w_accel *(xsi - oldxsi )
		xsi <- xsi + w_accel *(xsi - parm_history[,3] )
		parm_history[,1:2] <- parm_history[,2:3]
		parm_history[,3] <- xsi				
		xsi_change <- cbind( parm_history[,2] - parm_history[,1] , parm_history[,3] - parm_history[,2] )				
		lam <- eucl_norm( xsi_change[ind,2] )/ ( eucl_norm( xsi_change[ind,1] ) + eps )
		if ( iter > itermin ){
			w_accel <- lam / ( 2 - lam )
				if ( w_accel > 10 ){
						 w_accel <- 10
									}
								}
		xsi_acceleration$w <- w_accel  				
					}
	#*****************************************************
	# acceleration method of Ramsay (1975)
	if ( xsi_acceleration$acceleration == "Ramsay"){		  
       thre1 <- -5
	   # thre1 <- -2
		# use code from mirt package (author Phil Chalmers)				
        dX2 <- parm_history[ind,3] - parm_history[ind,2]
        dX <- xsi[ind] - parm_history[ind,3]
        d2X2 <- dX - dX2			
        # ratio <- sqrt((dX %*% dX) / (d2X2 %*% d2X2))
		ratio <- eucl_norm( dX ) / ( eucl_norm( d2X2 ) + eps )
        accel <- 1 - ratio
        if(accel < thre1){ accel <- thre1 }	
		if ( iter > itermin ){
			xsi <- (1 - accel) * xsi + accel * parm_history[,3]
							}				
		parm_history[,1:2] <- parm_history[,2:3]
		parm_history[,3] <- xsi			  		
		if ( iter > itermin ){		
			xsi_acceleration$beta_old <- accel				  
								}
							}
	#**************************************************************
	xsi_acceleration$parm_history <- parm_history
	xsi_acceleration$parm <- xsi	

	return(xsi_acceleration)	
						}
######################################################################						
# Euclidean norm of a vector
eucl_norm <- function(x){						
	sqrt( sum( x^2 ) )		
			}
######################################################################						
# Euclidean distance of two vectors			
eucl_distance <- function( x , y=NULL ){
	if ( is.matrix(x) ){
		y <- x[,2]
		x <- x[,1]
				}
	eucl_norm( x - y )
			}
######################################################################			