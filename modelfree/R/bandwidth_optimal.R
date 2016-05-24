bandwidth_optimal<-function( ptrue, r, m, x, H, link = c( "logit" ), guessing = 0,
                             lapsing = 0, K = 2, p = 1, ker = c( "dnorm" ),
                             maxiter = 50, tol = 1e-6, method = c( "all" ) ) {
#
# Finds the optimal bandwidth for a local polynomial estimate of the psychometric 
# function with specified guessing and lapsing rates. The difference between this function 
# and bandwidth_cross_validation is that here the true psychometric function is known.
# 
# INPUT
# 
# ptrue    - the true function; vector with the value of the true psychometric function at 
#		each stimulus level x
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
# K    - power parameter for Weibull and reverse Weibull link; default is 2
# p        - degree of the polynomial; default is 1
# ker      - kernel function for weights; default is "normpdf"
# maxiter  - maximum number of iterations in Fisher scoring; default is 50
# tol      - tolerance level at which to stop Fisher scoring; default is 1e-6
# method   - loss function to be used: choose from:
# 'ISEeta', 'ISE', 'deviance'; by default all possible values are calculated
#
# OUTPUT
# 
# h - optimal bandwidth for the chosen "method"; if no "method" is
# specified, then it has three components: $pscale, $eta-scale and $deviance 

# INTERNAL FUNCTIONS
# LOSS FUNCTION
    ISE <- function( f1, f2 ) {
        return( sum( ( f1 - f2 )^2 ) );
    }
# get ise for this value of h
# this is a nested function, so it shares variables!
    get_ise_p <- function( h ) {
        fest <- locglmfit( x, r, m, x, h, FALSE, link, guessing, lapsing,
                           K, p, ker, maxiter, tol );
# return MISE for this h
	return( ISE( ptrue, fest$pfit ) );
    }

# get ise on eta scale on for this value of h
# this is a nested function, so it shares variables!
    get_ise <- function( h ) {
        ftmp <- locglmfit( x, r, m, x, h, FALSE, link, guessing, lapsing,
                           K, p, ker, maxiter, tol );
        fit_eta <- linkfun( ptrue );
# return MISE for this h
        return( ISE( fit_eta, ftmp$eta ) );
    }
# get deviance for this value of h
# this is a nested function, so it shares variables!
    get_dev <- function( h ) {
        ftmp <- locglmfit( x, r, m, x, h, FALSE, link, guessing, lapsing,
                           K, p, ker, maxiter, tol ); 
# return MISE for this h
        D = return( deviance2( ptrue * m, m, ftmp$pfit ) );
    }

# MAIN PROGRAM
# First 5 arguments are mandatory
    if( missing("ptrue") || missing("r") || missing("m") || missing("x") ||
        missing("H") ) {
        stop("Check input. First 5 arguments are mandatory");
    }

# CHECK ROBUSTNESS OF INPUT PARAMETERS
# added error message in the case of length(ptrue) != length(x)
	if ( length(ptrue) != length(x)){
		stop("ptrue must have the same length as the stimulus levels x")
		}
    checkdata<-list();
    checkdata[[1]] <- x;
    checkdata[[2]] <- r;
    checkdata[[3]] <- m;
    checkinput( "psychometricdata", checkdata );
    rm( checkdata )
    checkinput( 'minmaxbandwidth', H );
    checkinput( 'linkfunction', link );
    if( length( guessing ) > 1 ) {
        stop( 'guessing rate must be scalar' );
    }
    if( length( lapsing ) > 1 ) {
        stop( 'lapsing rate must be scalar' );
    }
    checkinput( 'guessingandlapsing', c( guessing, lapsing ) );
    if (link == "weibull" || link == "revweibull"){
	    checkinput( "exponentk", K );
	    }
	pn <- list()
	pn[[1]] <- p
	pn[[2]] <- x
    checkinput( "degreepolynomial", pn );
    checkinput( 'kernel', ker );
    checkinput( 'maxiter', maxiter );
    checkinput( 'tolerance', tol );
    checkinput( 'method', method );

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

options(warn=-1)
# CROSS-VALIDATION
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
    options(warn=-1)
    return( h );
}
