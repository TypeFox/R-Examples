bandwidth_plugin<-function( r, m, x, link = c( "logit" ), guessing = 0,
                               lapsing = 0, K = 2, p = 1, ker = c( "dnorm" ) ) {
#
# The function calculates an estimate of the AMISE optimal bandwidth for 
# a local polynomial estimate of the psychometric function.
#
# INPUT
# 
# r    - number of successes at points x
# m    - number of trials at points x 
# x    - stimulus levels 
#
# OPTIONAL INPUT
# 
# link     - name of the link function to be used; default is "logit"
# guessing - guessing rate; default is 0
# lapsing  - lapsing rate; default is 0
# K    - power parameter for Weibull and reverse Weibull link; default is 2
# p        - degree of the polynomial; default is 1
# ker      - kernel function for weights; default is "dnorm"
#
# OUTPUT
# 
# h - plug-in bandwidth (ISE optimal on eta-scale)

# INTERNAL FUNCTIONS

# KERNELS
# Epanechnikov
    epanechnikov<-function( x ) {
        X <- x;
        X[which( abs( X ) > 1 )] <- 1;
        return( 0.75 * ( 1 - X^2 ) );
    }

# triangular
    triangular<-function( x ) {
        X <- x;
        X[which( abs( X ) > 1 )] <- 1;
        return( ( 1 - abs( X ) )  );
    }

# tri-cube
    tricube<-function( x ) {
        X <- x;
        X[which( abs( X ) > 1 )] <- 1;
        return( ( ( 1 - abs( X )^3 )^3 ) );
    }

# bi-square
    bisquare<-function( x ){
        X <- x;
        X[which( abs( X ) > 1 )] <- 1;
        return( ( (1 - abs( X )^2 )^2)  );
    }

# uniform
    uniform<-function( x ) {
        X <- x;
        return( ( abs( X ) <= 1 ) / 2 );
    }
    
    
# EQUIVALENT KERNEL
ker_eqv <- function( x, deg = p, K = ker, S1 = S ) {
    K2 <- NULL;
    K2 <- get( K );
    return( ( ( e_0 %*% solve( S1 ) %*% t( matrix( x, length( x ), deg + 1 )^
            t( matrix( 0:deg, deg + 1, length( x ) ) ) ) ) * K2( x ) ) );
}

ker_eqv2 <- function( x, deg = p, K = ker, S1 = S ) {
    return( ker_eqv( x, deg, K, S1 )^2 );
}

# MOMENTS OF KERNEL
    moments <- function( K, l ) {
        getxK <- function( x, power = 0, f = K ) {
             fun <- NULL;
             fun <- get( K );
             return( x^power * fun( x ) );
        }
        int <- integrate( getxK, -Inf, Inf, power = l );
        return( int$value );
    }

# VARIANCE FUNCTION FOR BINOMIAL DISTRIBUTION
    varfun <- function( x1, fit1 = fit ) {
# estimated mean
        mu <- as.numeric( predict( fit1, data.frame( x = x1 ), type = "response" ) );
# offset
        epsilon <- .001;
        ind <- rep( 0, length( x1 ) );
        ind[which( ( mu >= epsilon ) & ( mu <= 1 - epsilon ) )] <- 1;
# variance
        return( ind / ( mu  * (1 - mu ) ) );
    }

# MAIN PROGRAM
# First 3 arguments are mandatory
    if( missing("r") || missing("m") || missing("x") ) {
        stop("Check input. First 3 arguments are mandatory");
    }

# CHECK ROBUSTNESS OF INPUT PARAMETERS
    checkdata<-list();
    checkdata[[1]] <- x;
    checkdata[[2]] <- r;
    checkdata[[3]] <- m;
    checkinput( "psychometricdata", checkdata );
    rm( checkdata )
    pn <- list()
	pn[[1]] <- p
	pn[[2]] <- x
    checkinput( "degreepolynomial", pn );
    checkinput( "kernel", ker );
    checkinput( 'linkfunction', link );
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
    
    n <- length(r);

# p+Dp - degree of polynomial used for parametric fit
    if( ( n - p ) > 5 ) {
        Dp <- 3;
    }
    else {
        if( ( n - p ) == 5 ) {
            Dp <- 2;
        }
        else {
            if ( ( n - p ) > 2 ) {
                Dp = 1; 
            }
            else {
                stop( paste( "Not enough data to fit polynomial of degree", p,
                      sep = " " ) );
            }
        }
    }

# PARAMETRIC ESTIMATOR
# parametric estimator of order p + Dp

# INITIAL
    p1 <- p + Dp;
    lx <- length( x );

# GLOBAL FIT
# create data frame
    glmdata <- data.frame( cbind( r/m ,m , x ) );
    names( glmdata ) <- c( "resp", "m", "x" );
# formula
    glmformula <- c( "resp ~ x" );
    if( p1 > 1 ) {
        for( pp in 2:p1 ) {
            glmformula <- paste( glmformula, " + I(x^", pp,")", sep = "");
        }
    }
# fit
    if( ( link == "weibull"    ) || 
        ( link == "revweibull" ) ) {
        linkfun <- get( paste( link, "_link_private", sep = "" ) );
        fit <- glm( glmformula, data = glmdata, weights = m,
                    family = binomial( linkfun( K, guessing, lapsing ) ) );
    }
    else {
    	if( link == "logit"      ||
        link == "probit"     ||
        link == "loglog"     ||
        link == "comploglog" ){
        	
        	linkfun <- get( paste( link, "_link_private", sep = "" ) );
        	fit <- glm( glmformula, data = glmdata, weights = m,
                    family = binomial( linkfun( guessing, lapsing ) ) );
                    }
                    else{
                    	linkfun <- get( link );
        				fit <- glm( glmformula, data = glmdata, weights = m,
                    			family = binomial( linkfun( guessing, lapsing ) ) );
                    	}
    }


if (any(!is.finite(fit$coefficients)) ) stop('Result is degenerate: the link function is probably inappropriate')
    pol <- polynom( fit$coefficients );

# adjust degenrated values to aviod degenerate results
    epsilon <- .001;
    ind <- rep( 0, lx );
    ind[which( ( r / m >= epsilon ) & ( r / m <= 1 - epsilon ) )] <- 1;

# coefficients for (p+1)th derivative of parametric estimator
    for( l in 1:(p+1) ) {
        pol <- deriv( pol );
    }
    tmp_b_p1 <- coefficients( pol );

# calculate the (p+1)th derivative of eta
    tmp_eta_p1 <-pol( x ) * ind;

# approximate int(eta^(p+1))^2
    int_eta <- sum( tmp_eta_p1^2 ) * m[1];

# EQUIVALENT KERNEL
    tmp <- NULL;
    for( s in 0:(2 * p) ) {
        tmp[s+1] <- moments( ker, s );
    }

# create S matrix (all moments)
    S <- matrix( 0, p + 1, p + 1 );
    for(pp in 1:(p+1) ) {
        S[pp,] <- tmp[pp+(0:p)];
    }
# indicator vector
    e_0 <- diag( 1, 1, p + 1 );

# FUNCTIONS OF THE EQUIVALENT KERNEL
# integral of squared equivalent kernel
    k2 <- integrate( ker_eqv2, -Inf, Inf, deg = p, K = ker )$value;
# (p+1)th moment of equivalent kernel 
    muK_2 <- moments( "ker_eqv", p + 1 )^2;
    
 if( muK_2 == 0 ) {
        muK_2 <- .Machine$double.eps;
    }
# INTEGRAL OF VARIANCE FUNCTION
    int_sg <- integrate( varfun, min( x ) - 1, max( x ) + 1 )$value;
# ensure the above integral is positive which might not be the case for
# near degenerate samples; issue warning 
    if( int_sg == 0 ) {
        int_sg <- .Machine$double.eps;
        warning("The estimated value was 0; sample is probably degenerate");
    }
# KERNEL DEPENDANT CONSTANT
    C <- ( ( factorial( p + 1 )^2 * k2 ) / ( 2 * ( p + 1 ) * muK_2 ) )^
         ( 1 / ( 2 * p + 3 ) );


# BANDWIDTH
    return( C * ( int_sg / int_eta )^( 1 / ( 2 * p + 3 ) ) );
}
