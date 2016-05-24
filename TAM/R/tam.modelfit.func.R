
########################################################
# tam.modelfit with user defined input
tam.modelfit.args <- function( resp , probs , theta , post , progress=TRUE ){	    
		resp.ind <- as.matrix( 1- is.na(resp) )
		tamobj <- list( "resp" = resp , "rprobs" = probs , 
					"theta" = theta , "hwt" = post ,
					"resp.ind" = resp.ind )				
		res <- tam.modelfit( tamobj=tamobj , progress=progress)			
		return(res)																						
		}
########################################################		

###############################################################	
tam.modelfit.IRT <- function( object , progress=TRUE ){	
	resp <- IRT.data(object)
	probs <- IRT.irfprob(object)
	theta <- attr( probs , "theta" )
	post <- IRT.posterior( object )
	res <- tam.modelfit.args( resp , probs , theta , post , progress)
	return(res)
		}
#############################################################