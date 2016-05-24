
###################################
#*** delta design matrix
tamaanify.tam.mml.3pl.deltadesign <- function(res){
	anlist <- res$ANALYSIS.list
	delta.fixed <- NULL
	#*******************************************
	#*** delta design matrix mixture distribution
	if ( res$ANALYSIS.list$type == "MIXTURE" ){
		theta <- res$theta_MIXTURE
		TG <- nrow(theta)
		ncl <- anlist$NCLASSES
		D <- ncol(theta)
		TP <- TG*ncl
		delta.designmatrix <- matrix( 0 , nrow=TP , ncol=0)
		th1 <- create.deltaDesign(X=theta , main=2 , int=TRUE )
		# bind two matrices
		delta.designmatrix <- th1
		for (cl in 2:ncl){
			delta.designmatrix <- rbind.twomatrices(X1=delta.designmatrix , X2=th1)
						}
		rownames(delta.designmatrix) <- dimnames(res$E)[[3]]
		
# Revalpr("delta.designmatrix")
# delta.designmatrix <- NULL		
		
	    res$delta.designmatrix <- delta.designmatrix
		
#		delta.fixed <- cbind( c(2,5) , 1 , 0 )
		
					}
	#******************************************
	
	res$delta.fixed <- delta.fixed
	
	return(res)
	}