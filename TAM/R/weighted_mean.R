############################################
# weighted mean
# wrapper to weighted.mean
weighted_mean <- function( x , w=rep(1,length(x)) ){
	res <- stats::weighted.mean(x=x , w=w , na.rm=TRUE )
	return(res)
		}
###############################################