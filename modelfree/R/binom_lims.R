binom_lims<-function( r, m, x, gl = c( "both" ), link = c( "logit" ), p = 1, K = 2,
                      initval = NULL ) {
#
# This function finds the maximum likelihood estimates of the parameters 
# of the psychometric function with guessing and lapsing rates, only 
# guessing rate, or only lapsing rate.  
#
# INPUT
#
# r    - number of successes at points x
# m    - number of trials at points x 
# x    - stimulus levels
#
# OPTIONAL INPUT
# 
# gl      - indicator, calulate only guessing if "guessing", only lapsing if "lapsing"
#  and both guessing and lapsing if "both"; default is "both"
# link    - name of the link function; default is "logit"
# p       - degree of the polynomial; default is 1
# K       - power parameter for Weibull and reverse Weibull link; default is 2
# initval - initial value for guessing and lapsing; default is c(.01 .01) if guessing and
# lapsing rates are estimated, and .01 if only guessing or only lapsing rate is estimated
#
# OUTPUT
# 
# Object with 3 or 4 components: 
# b - estiamted coefficients for the linear part
# guessing - estimated guessing rate (if estimated)
# lapsing - estimated lapsing rate (if estimated)
# fit - glm object to be used in evaluation of fitted values

# MAIN PROGRAM
# First 3 arguments are mandatory
    if( missing("x") || missing("r") || missing("m") ) {
        stop("Check input. First 3 arguments are mandatory");
    }

    if( gl != "both"  && gl != "guessing" && gl != "lapsing") {
        stop( "Wrong value for guessing/lapsing indicator gl" );
    }
    if(is.null(initval)){
    
    if( gl == "guessing" || gl == "lapsing" ){ 
    	initval <- .01;
    	}else{
    		initval <- c(.01,.01);
    		} 
    }

# CHECK ROBUSTNESS OF INPUT PARAMETERS
    checkdata<-list();
    checkdata[[1]] <- x;
    checkdata[[2]] <- r;
    checkdata[[3]] <- m;
    checkinput( "psychometricdata", checkdata );
    rm( checkdata )
    checkinput( "linkfunction", link );
    pn <- list()
	pn[[1]] <- p
	pn[[2]] <- x
    checkinput( "degreepolynomial", pn );
    if( link == "weibull"  || link == "revweibull") {
    	checkinput( 'exponentk', K );
    	}
    	
    initval[initval==0] <- .Machine$double.eps
    
   
# Check initial values for guessing or guessing/lapsing and call internal
# functions binom_gl or binom_g depending on indicator gl
    if( gl == "both" ) {
    	if( length(initval) != 2){
    		stop( "initval must have two elements if both guessing and lapsing rates are being estimated")
    		}
        checkinput( "guessingandlapsing", initval );
        return( binom_gl( r, m, x, link, p, K, initval ) );
    }
    else {
    	if( gl == "guessing" ) {
# Check that initval is a positive scalar
        if( length( initval ) > 1 || initval <= 0 || initval >= 1 ) {
            stop( "Guessing rate must be a scalar between 0 and 1" );
        }
        return( binom_g( r, m, x, link, p, K, initval ) );
        }# end guessing
        else{
        	if( length( initval ) > 1 || initval <= 0 || initval >= 1 ) {
            stop( "Lapsing rate must be a scalar between 0 and 1" );
	        }
    	    return( binom_l( r, m, x, link, p, K, initval ) );
        	}# end lapsing

        }
    
}
