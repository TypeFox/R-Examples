
#########################################################################
# computation information criteria
rm.facets.IC <- function( dev , dat2 , VV , RR , maxK , a.item.center,
		est.a.item , est.b.rater , est.a.rater , est.mean ,
		b.rater.center , a.rater.center , b.rater.fixed , 
		a.rater.fixed , tau.item.fixed_val , a.item.fixed
		){
	# Information criteria
	ic <- list( "deviance" = dev , "n" = nrow(dat2) )
	ic$VV <- VV
	ic$RR <- RR
	# item parameters
	# tau.item
	ic$np.item <- sum(maxK)
	if ( ! is.null( tau.item.fixed_val ) ){
	   ic$np.item <- ic$np.item - sum( 1-is.na( tau.item.fixed_val ) )
										}
	# a.item
	b1 <- VV-a.item.center
	if ( ! is.null(a.item.fixed) ){
	    b1 <- b1 - sum( ! is.na( a.item.fixed ) )
									}	
	ic$np.item <- ic$np.item + est.a.item*b1
	#****
	#-- rater parameters
	ic$np.rater <- 0	
	# b.rater
	b1 <- RR-b.rater.center
	if ( ! is.null(b.rater.fixed) ){
	    b1 <- b1 - sum( ! is.na( b.rater.fixed ) )
									}
	ic$np.rater <- est.b.rater*b1
    # a.rater	
	b1 <- RR-a.rater.center
	if ( ! is.null(a.rater.fixed) ){
	    b1 <- b1 - sum( ! is.na( a.rater.fixed ) )
									}	
	ic$np.rater <- ic$np.rater + est.a.rater*b1
	# distribution parameters
	ic$np.trait <- 1 + est.mean
	# estimated parameters
	ic$np <- ic$np.trait + ic$np.item + ic$np.rater
    # AIC
    ic$AIC <- dev + 2*ic$np
    # BIC
    ic$BIC <- dev + ( log(ic$n) )*ic$np
    # CAIC (consistent AIC)
    ic$CAIC <- dev + ( log(ic$n) + 1 )*ic$np
	# corrected AIC
    ic$AICc <- ic$AIC + 2*ic$np * ( ic$np + 1 ) / ( ic$n - ic$np - 1 )	
	return(ic)
	  }
#########################################################################	