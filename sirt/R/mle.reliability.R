

#########################################################
# MLE reliability
mle.reliability <- function(meas, se.meas ){
		meas[ abs(meas) > 1E10 ] <- NA
		v1 <- stats::var(meas , na.rm=TRUE)
		v2 <- mean( se.meas^2 , na.rm=TRUE )		
		rel <- ( v1 - v2 ) / v1
		return(rel)
			}
#########################################################			

#########################################################
# EAP reliability
eap.reliability <- function(meas, se.meas ){
		meas[ abs(meas) > 1E10 ] <- NA
		v1 <- stats::var(meas , na.rm=TRUE)
		v2 <- mean( se.meas^2 , na.rm=TRUE )		
		rel <- v1  / ( v1 + v2 )
		return(rel)
			}
#########################################################		