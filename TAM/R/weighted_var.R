
#####################################################################
# weighted variance.
# This function is just a wrapper to cov.wt
weighted_var <- function( x , w=rep(1,length(x) ) , method = "unbiased"  ){
    ind <- which( ! is.na(x) )
    x <- x[ind]
    w <- w[ind]
	dat <- data.frame("x" = x )		
	res <- stats::cov.wt( x = dat , wt = w , cor =FALSE , center = TRUE ,
	            method = method )
	res <- res$cov[1,1]
	return(res)	
	}
#######################################################################	
# standard deviation
weighted_sd <- function( x , w=rep(1,length(x) ) , method = "unbiased"  ){
	res <- sqrt( weighted_var(x=x, w=w, method=method) )
	return(res)
			}
#######################################################################