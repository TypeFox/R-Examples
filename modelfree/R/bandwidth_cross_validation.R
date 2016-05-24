bandwidth_cross_validation<-function( r, m, x, H, link = c( "logit" ),
                                      guessing = 0, lapsing = 0, K = 2,
                                      p = 1, ker = c( "dnorm" ), maxiter = 50,
                                      tol = 1e-6, method = c( "all" ) ) {
#
# The function finds the cross-validation bandwidth for a local polynomial estimate of 
# the psychometric function with specified guessing and lapsing rates.
#
# INPUT
# 
# r    - number of successes at points x
# m    - number of trials at points x 
# x    - stimulus levels 
# H    - search interval
#
# OPTIONAL INPUT
#
# link     - name of the link function to be used; default is "logit"
# guessing - guessing rate; default is 0
# lapsing  - lapsing rate; default is 0
# K        - power parameter for Weibull and reverse Weibull link; default is 2
# p        - degree of the polynomial; default is 1
# ker      - kernel function for weights; default is "dnorm"
# maxiter  - maximum number of iterations in Fisher scoring; default is 50
# tol      - tolerance level at which to stop Fisher scoring; default is 1e-6
# method   - loss function to be used in cross-validation: choose from:
# "ISEeta", "ISE", "deviance"; by default all possible values are calculated
#
# OUTPUT
# 
# h - cross-validation bandwidth for the chosen "method"; if no "method" is
# specified, then it has three components: $pscale, $etascale and $deviance 

#### 
# KZ 28-Mar-12
# included on.exit function which restores warning settings to their
# original state
####


# INTERNAL FUNCTIONS
# LOSS FUNCTION
    ISE <- function( f1, f2 ) {
        return( sum( ( f1 - f2 )^2 ) );
    }
# get ise for this value of h
# this is a nested function, so it shares variables!
    get_ise_p <- function( h ) {
        fest <- NULL;
        for( i in 1:Lx ) {
 
 		    fest[i] <- locglmfit( x[i], r[-i], m[-i], x[-i], h, FALSE,
                        link, guessing, lapsing, K, p, ker, maxiter, tol )$
                        pfit;
                    
        }
 
 # return MISE for this h

	return( ISE( r / m, fest ) );
    }

# get ise on eta scale on for this value of h
# this is a nested function, so it shares variables!
    get_ise <- function( h ) {

        fest <- NULL;
        for( i in 1:Lx ) {
            fest[i] <- locglmfit( x[i], r[-i], m[-i], x[-i], h, FALSE,
                       link, guessing, lapsing, K, p, ker, maxiter, tol )$
                       etafit;
        }
        fit <- ( r + .5 ) / ( m + 1 );
        fit_eta <- linkfun( fit )
        
# return MISE for this h
        return( ISE( fit_eta, fest ) );
    }
# get deviance for this value of h
# this is a nested function, so it shares variables!
    get_dev <- function( h ) {

        fest <- NULL;
        for( i in 1:Lx ) {

            fest[i] <- locglmfit( x[i], r[-i], m[-i], x[-i], h, FALSE,
                       link, guessing, lapsing, K, p, ker, maxiter, tol )$
                       pfit;
        }
# return MISE for this h
        return( deviance2( r, m, fest ) );
    }

# MAIN PROGRAM
# First 4 arguments are mandatory
    if( missing("r") || missing("m") || missing("x") || missing("H") ) {
        stop("Check input. First 4 arguments are mandatory");
    }

# CHECK ROBUSTNESS OF INPUT PARAMETERS

    checkdata<-list();
    checkdata[[1]] <- x;
    checkdata[[2]] <- r;
    checkdata[[3]] <- m;
    checkinput( "psychometricdata", checkdata );
    rm( checkdata )
    checkinput( "minmaxbandwidth", H );
    checkinput( "linkfunction", link );
    if( length( guessing ) > 1 ) {
        stop( "guessing rate must be scalar" );
    }
    if( length( lapsing ) > 1 ) {
        stop( "lapsing rate must be scalar" );
    }
    checkinput( "guessingandlapsing", c( guessing, lapsing ) );
    checkinput( "guessingandlapsing", c( guessing, lapsing ) );
    if (link == "weibull" || link == "revweibull"){
	    checkinput( "exponentk", K );
	    }
	pn <- list()
	pn[[1]] <- p
	pn[[2]] <- x
    checkinput( "degreepolynomial", pn );
    checkinput( "kernel", ker );
    checkinput( "maxiter", maxiter );
    checkinput( "tolerance", tol );
    checkinput( "method", method );# initial values for the loops
    LH = length(H);
    Lx = length(x);
    
    # Retrieve link and inverse link functions
    if( link == "logit"      ||
        link == "probit"     ||
        link == "loglog"     ||
        link == "comploglog" ||
        link == "weibull"    ||
        link == "revweibull" ) {
            	
    			LINK <- paste( link, "_link_private", sep = "" );
	}else{
		LINK <- link
		}

    if( LINK != "weibull_link_private" && LINK != "revweibull_link_private" ) {

        linkuser <- eval( call( LINK, guessing, lapsing ) );
        
    }
    else{

       	linkuser <- eval( call( LINK, K,  guessing, lapsing ) );

    }
      
    linkfun <- linkuser$linkfun;


# KZ 28-03-2012 included on.exit routine so that the warning settings are
# restored when the function terminates even if interrupted by user

warn.current <- getOption("warn")
on.exit(options(warn = warn.current));

options(warn=-1)


# BANDIWDTH
    h <- NULL;
    if( method == "ISE" ) {
# p-scale
        h <- optimize( get_ise_p, lower = H[1], upper = H[2] )$minimum;
     }
     else {
         if ( method == "ISEeta") {
# eta scale
             h <- optimize( get_ise, lower = H[1], upper = H[2] )$minimum;
         }
         else {
             if( method == "likelihood" ) {
# DEVIANCE
                 h <- optimize( get_dev, lower = H[1], upper = H[2] )$minimum;
             }
             else {
# p-scale
                 h$pscale   <- optimize( get_ise_p, lower = H[1], upper = H[2]
                               )$minimum;
# eta scale
                 h$etascale <- optimize( get_ise, lower = H[1], upper = H[2]
                               )$minimum;

# DEVIANCE
                 h$deviance <- optimize( get_dev, lower = H[1], upper = H[2]
                               )$minimum;
            }
        }
    }
	

    return( h );
}
