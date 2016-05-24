bandwidth_bootstrap <- function( r, m, x, H, N, h0 = NULL, link = c("logit"), guessing = 0,
                               lapsing = 0, K = 2, p = 1,
                               ker = c("dnorm"), maxiter = 50, tol = 1e-6,
                               method = c("all") ) {
#
# The function finds a bootstrap estimate of the optimal bandwidth h for a local polynomial 
# estimate of the psychometric function with specified guessing and lapsing rates.
#
# INPUT
# 
# r    - number of successes at points x
# m    - number of trials at points x 
# x    - stimulus levels 
# H    - search interval
# N    - number of bootstrap replications
#
# OPTIONAL INPUT
#
# h0       - pilot bandwidth; if not specified, then the scaled plug-in
# 			bandwidth is used
# link     - name of the link function to be used; default is "logit"
# guessing - guessing rate; default is 0
# lapsing  - lapsing rate; default is 0
# K    - power parameter for Weibull and reverse Weibull link; default is 2
# p        - degree of the polynomial; default is 1
# ker      - kernel function for weights; default is "dnorm"
# maxiter  - maximum number of iterations in Fisher scoring; default is 50
# tol      - tolerance level at which to stop Fisher scoring; default is 1e-6
# method   - loss function to be used in bootstrap: choose from:
# "ISEeta", "ISE", "deviance"; by default all possible values are calculated
#
# OUTPUT
# 
# h - bootstrap bandwidth for the chosen "method"; if no "method" is
# specified, then it has three components: $pscale, $eta-scale and $deviance

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
        fest <- matrix( 0, n, N );
        for( i in 1:N ) {
            fest[,i] <- locglmfit( x, samp[,i], m, x, h, FALSE, link,
                        guessing, lapsing, K, p, ker, maxiter, tol )$pfit;

}
# return MISE for this h
		return( mean( ISE( f1, fest ) ) );
    }
# get ise on eta scale on for this value of h
# this is a nested function, so it shares variables!
    get_ise <- function( h ) {
        fest <- matrix( 0, n, N );
        for( i in 1:N ) {
            fest[,i] <- locglmfit( x, samp[,i], m, x, h, FALSE, link,
                        guessing, lapsing, K, p, ker, maxiter, tol )$etafit;
        }
# return MISE for this h
	return( mean( ISE( eta1, fest ) ) );
    }
# get deviance for this value of h
# this is a nested function, so it shares variables!
    get_dev <- function( h ) {
        devfit <- numeric( N );
        for( i in 1:N ) {
            fest <- locglmfit( x, samp[,i], m, x, h, FALSE, link,
                    guessing, lapsing, K, p, ker, maxiter, tol )$pfit;
            devfit[i] <- deviance2( r, m , fest );
        }
# return MISE for this h
        return( mean( devfit ) );
    }

# MAIN PROGRAM
# First 5 arguments are mandatory
    if( missing("r") || missing("m") || missing("x") || missing("H") ||
        missing("N")) {
        stop("Check input. First 5 arguments are mandatory");
    }
# CHECK ROBUSTNESS OF INPUT PARAMETERS
    checkdata<-list();
    checkdata[[1]] <- x;
    checkdata[[2]] <- r;
    checkdata[[3]] <- m;
    checkinput( "psychometricdata", checkdata );
    rm( checkdata )
    checkinput( "minmaxbandwidth", H );
    checkinput( "bootstrapreplications", N );
    if( !is.null( h0 ) ) {
        checkinput( "bandwidth", h0 );
    }
    
    checkinput( "linkfunction", link );
    if( length( guessing ) > 1 ) {
        stop( "Guessing rate must be a scalar" );
    }
    if( length( lapsing ) > 1 ) {
        stop( "Lapsing rate must be a scalar" );
    }
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
    checkinput( "method", method );

    n <- length(x);

# KZ 28-03-2012 included on.exit routine so that the warning settings are
# restored when the function terminates even if interrupted by user

warn.current <- getOption("warn")
on.exit(options(warn = warn.current));

options(warn=-1)


# OBTAIN INITIAL BANDWIDTH, IF NOT GIVEN
    if( is.null( h0 ) ) {
        h0 <- (1.5 * n^.1) * bandwidth_plugin( r, m, x, link, guessing, lapsing, K, p, ker );
    }

# INITIAL ESTIMATE
# pilot over-smoothed estimate with bandiwdth h0
    f <- locglmfit( x, r, m, x, h0, FALSE, link, guessing, lapsing, K,
                    p, ker, maxiter, tol );

# SAMPLING
# re-sampling
    samp <- matrix( 0, n, N );
	samp <- matrix(rbinom(N*n,m,f$pfit),n,N)

# exclude "degenerate samples" if min(M)>1
    if( min( m ) > 1 ) {
        ok <- NULL;
        for( i in 1:N ) {
            ok[i] = length( unique( samp[,i] ) ) > 3;
        }
        while( any( ok == FALSE ) ) {
            lis <- which( ok == 0 );
  			samp[,lis] <- matrix(rbinom(length(lis)*n,m,f$pfit),n,length(lis))
            for( i in 1:length( lis ) ) {
                ok[lis[i]] <- length( unique( samp[,lis[i]] ) ) > 3;
            }
        }
    }

# INITIATE VARIABLE IN WHICH DATA ARE STORED
f1 <- matrix( rep( f$pfit, N ), n, N );
eta1 <- matrix( rep( f$etafit, N ), n, N );

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