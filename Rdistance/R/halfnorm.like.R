halfnorm.like <- function(a, dist, w.lo=0, w.hi=max(dist), series="cosine", expansions=0, scale=TRUE){
#   Computes half norm likelihood, scaled appropriately to integrate to 1.0, for every
#   observation. I.e., returns a vector. 
#
#   Input:
#   a = parameter values. Length and meaning depend on series and expansions
#   dist = input observed distance data
#   w = right truncation value, same units as dist
#   series = character values specifying type of expansion.  Currently only 
#       "cosine" and "hermite" work. Default is no series expansions.
#   expansions = number of expansion terms.  This controls whether series 
#       expansion terms are fitted. Default = 0 does not fit any terms. 
#       If >0, fit the number of expansion terms specified of the type 
#       specified by series.  Max terms depends on series.
#   scale = logical, whether or not to scale the likelihood values so 
#       that the underlying density integrates to 1. This parameter is used 
#       to stop the recursion. 
#   
#   Output:
#   A vector same length as dist, containing likelihood values, scaled appropriately. 
#   Likelihood for all distances >w are set to NA
#


    dist[ (dist < w.lo) | (dist > w.hi) ] <- NA

    sigma <- a[1]
    key <- exp(-dist^2/(2*sigma^2))
    dfunc <- key
    w <- w.hi - w.lo
	
    # If there are expansion terms
    if(expansions > 0){

        nexp <- min(expansions,length(a)-1)  # should be equal. If not, fire warning next
        
        if( length(a) != (expansions+1) ) {
            warning("Wrong number of parameters in expansion. Should be (expansions+1). High terms ignored.")
        }

		if (series=="cosine"){
            dscl = dist/w
            exp.term <- cosine.expansion( dscl, nexp )
		} else if (series=="hermite"){
            dscl = dist/sigma
            exp.term <- hermite.expansion( dscl, nexp )
		} else if (series == "simple") {
            dscl = dist/w
            exp.term <- simple.expansion( dscl, nexp )
        } else {
            stop( paste( "Unknown expansion series", series ))
        }

        dfunc <- key * (1 + (exp.term %*% a[2:(nexp+1)]))


    } else if(length(a) > 1){
        warning("Wrong number of parameters in halfnorm. Only 1 needed if no expansions. High terms ignored.")
    }

    if( scale ){
        dfunc = dfunc / integration.constant(halfnorm.like,w.lo=w.lo,w.hi=w.hi,a=a,series=series,expansions=expansions)   # scales underlying density to integrate to 1.0
        
           #df2 <- dfunc[ order(dist) ]
           #d2 <- dist[ order(dist) ]
            #cat(paste("integral=", sum( diff(d2) * (df2[-1] + df2[-length(df2)]) ) / 2, "\n" ))
        
            #readline("Enter:")

    }
    
    

    c(dfunc)
    
}
