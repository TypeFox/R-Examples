negexp.like <- function (a, dist, w.lo=0, w.hi=max(dist), series="cosine", expansions = 0, scale=TRUE){
#
#   Compute negative exponential likelihood
#
#   Input:
#   a = parameter values. Length and meaning depend on series and expansions.
#   dist = input observed distance data
#   w = right truncation value, same units as dist
#   series = character values specifying type of expansion.  Currently only
#       "cosine" and "simple" work. Default is no series expansions.
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
#   Details: 	
#   Excerpted from DISTANCE online help manual:

#   The candidate key functions offered are Uniform, Half-normal, Hazard-rate and
#   Negative exponential (this last function is not recommended, except for salvage
#   analyses; see Buckland et al. 1993, 2001). The candidate series expansions are
#   Cosine, Simple polynomial and Hermite polynomial.
#
#   Bounds on Key Function Parameters
#   With some datasets, it may be necessary to bound the parameter estimates to
#   achieve convergence. In these situations, you can use this table to specify upper
#   and lower bounds on the key function parameters
#
#   (one parameter for half-normal and negative exponential key functions,
##--THIS IS WHERE I GOT THE ONE PARAMETER NEG EXP
#
#   two for the hazard rate function, plus any covariate parameters). One common
#   circumstance where this is required is to impose a lower bound of 1.0 on the
#   second hazard rate parameter

    dist[ (dist < w.lo) | (dist > w.hi) ] <- NA

	beta=a[1]
	key = exp(-beta*dist)
    dfunc <- key
    w <- w.hi - w.lo


    if(expansions > 0){

        nexp <- min(expansions,length(a)-1)  # should be equal. If not, fire warning next
        
        if( length(a) != (expansions+1) ) {
            warning("Wrong number of parameters in expansion. Should be (expansions+1). High terms ignored.")
        }

		if (series=="cosine"){
            dscl = dist/w
            exp.term <- cosine.expansion( dscl, nexp )
		} else if (series=="hermite"){
            dscl = dist/w
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
        dfunc = dfunc / integration.constant(negexp.like,w.lo=w.lo,w.hi=w.hi,a=a,series=series,expansions=expansions)  # makes integral from w.lo to w.hi = 1.0
    }
    
    c(dfunc)
}
