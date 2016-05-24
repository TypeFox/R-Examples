
#####################################################
# skewness
weighted_skewness <- function( x , w = rep(1,length(x)) ){	
	m <- weighted_mean( x=x , w=w)
	s <- weighted_sd( x=x, w=w)
	y <- ((x-m)/s)^3
	res <- weighted_mean(x=y, w=w)	
	return(res)
		}
#######################################################		
# curtosis
weighted_curtosis <- function( x , w = rep(1,length(x)) ){	
	m <- weighted_mean( x=x , w=w)
	v <- weighted_var( x=x, w=w)
	y <- (x-m)^4
	res <- weighted_mean(x=y, w=w) / v^2 - 3	
	return(res)
		}
###########################################################		