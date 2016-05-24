#####
## This function from risksetROC
#####
IntegrateAUC <- function( AUC, utimes, St, tmax, weight="rescale" )
{
##
## Assume ordered utimes
##
	wChoice <- match( weight, c("rescale","conditional") )
	
	if( is.na(wChoice) )
    {
		cat("error in weight choice")
		stop(0)
	}
	ft <- rep( NA, length(St) )
	ft[1] <- 1.0 - St[1]
## this is original line:
## for( j in 1:length(St) ) ft[j] <- St[j-1] - St[j]
## this is corrected line: (?)
	for( j in 2:length(St) ) ft[j] <- St[j-1] - St[j]
##
	mIndex <- length( utimes[ utimes <= tmax ] )
	www <- 2*ft*St
	iAUC <-  0.0
	wTotal <- sum( www[1:mIndex] )
	Smax <- St[ min(mIndex+1,length(St)) ]
##   for( j in 1:mIndex ){
##     if( wChoice==1 ){
##       wj <- 2*ft[j]*St[j]/wTotal
##     }
##     if( wChoice==2 ){
##       wj <- 2*( ft[j]/(1-Smax) )*(St[j]-Smax)/(1-Smax)
##     }
##     iAUC <- iAUC + wj*AUC[j]
##   }
##   iAUC
##   
	if(wChoice == 1)      w=2*ft*St/wTotal
	else  w = 2*ft*(St-Smax)/((1-Smax)^2)
	iAUC=sum(w[1:mIndex]*AUC[1:mIndex])
	return(iAUC)
}