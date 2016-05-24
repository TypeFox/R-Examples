
#################################
# algorithm for Aicken's statistic
agree_aicken <- function( PAk , PBk , Pa ){
	PAH <- PAk
	PBH <- PBk
	alpha <- .5
	ii <- 0
	conv <- FALSE
	maxiter <- 100
	globconv <- 1E-5
	while( ! conv){
		alpha0 <- alpha
		pet <- sum( PAH * PBH )
		alpha <- ( Pa - pet ) / ( 1 - pet )
		PAH <- PAk / (  ( 1 - alpha )  + alpha * PBH  / pet )
		PBH <- PBk / (  ( 1 - alpha )  + alpha * PAH  / pet )
		ii <- ii+1
		diff_conv <- abs( alpha0 - alpha )
		if (diff_conv < globconv){ conv <- TRUE }
		if (ii == maxiter ){ conv <- TRUE }
			}
	# chance agreement
	Pe <- pet
	# output
	res <- list( "alpha" = alpha , "PAH" = PAH ,
			"PBH" = PBH , "Pe" = Pe )
	return(res)
		}