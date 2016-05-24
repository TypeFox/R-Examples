
###################################################
# summary for tam.fit
summary.tam.fit <- function( object , ... ){
    object <- object$itemfit
	ind <- grep( "pholm" , colnames(object) )
	obji <- object[ , - ind ]
	for ( vv in seq(2,ncol(obji) ) ){
		obji[,vv] <- round( obji[,vv] , 3 )
					}
	return(obji)
		}
###################################################