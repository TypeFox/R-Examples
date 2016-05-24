
#####################################################
# computes reliability for one-dimensional WLEs
WLErel <- function( theta , error , w = rep(1,length(theta) )){
    pweights <- w		
    v1 <- weighted_var( x = theta , w = pweights  )
    v2 <- stats::weighted.mean( error^2 , w=pweights , na.rm=TRUE )	
    # WLE_Rel = ( v1 - v2 ) / v1 = 1 - v2 / v1
    WLE.rel <- 1 - v2 / v1
	return(WLE.rel)
	}
#######################################################	
EAPrel <- function( theta , error , w = rep(1,length(theta) )){
    pweights <- w		
    v1 <- weighted_var( x = theta , w = pweights  )
    v2 <- stats::weighted.mean( error^2 , w=pweights , na.rm=TRUE )	
    # v1 / (v1+v2) = 1 - v2 / ( v1 + v2 )
    rel <- v1 / ( v1 + v2 )
	return(rel)
	}
#######################################################	
				