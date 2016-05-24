binom_l<-function( r, m, x, link, p, K,
                      initval ) {
#
# THIS IS AN INTERNAL FUNCTION: USE BINOM_LIMS FOR BEST RESULTS
#
# Maximum likelihood estimates of the parameters of the psychometric 
# function with lapsing rate. The estimated parameters for the linear 
# part are in vector 'b' and the estimated lapsing rate is 'lapsing'.
#
# INPUT
#
# r    - number of successes at points x
# m    - number of trials at points x 
# x    - stimulus levels
# link    - name of the link function
# p       - degree of the polynomial
# K       - power parameter for Weibull and reverse Weibull link
# initval - initial value for lapsing
#
# OUTPUT
# 
# Object with 3 components: 
# b - vector of estiamted coefficients for the linear part 
# lapsing - estimated lapsing rate
# fit - glm object to be used in evaluation of fitted values

# LIKELIHOOD
    likfun <- function( lapsing ) {

        lapsing <- 1 / ( 1 + exp( -lapsing ) );
        guessing <- 0;

# fit
        if( linkfun != "weibull_link_private" && linkfun != "revweibull_link_private" ) {
            fit <- glm( glmformula, data = glmdata, weights = m,
                        family = binomial( eval( call( linkfun, guessing, lapsing ) ) ) );

        }
        else {
            fit <- glm( glmformula, data = glmdata, weights = m,
                        family = binomial( eval( call( linkfun, K, guessing, lapsing ) ) ) );
        }
# FITTED PROBABILITIES
        fitted <- as.numeric( predict( fit, type = "response" ) );
        fitted[which( fitted <= 0 )] <- 0 + .Machine$double.eps;
        fitted[which( fitted >= 1-lapsing )] <- 1-lapsing - .Machine$double.eps;

        return( c( -( t( r ) %*% log( fitted ) + t( m - r ) %*%
                log( 1 - fitted ) ) ) );
                
    }

# MAIN PROGRAM

    initval <- log( initval / ( 1 - initval ) );
    
# GLM settings
    glmdata <- data.frame( cbind( r/m , m , x ) );
    names( glmdata ) <- c( "resp", "m", "x" );

# formula
    glmformula <- c( "resp ~ x" );
    if( p > 1 ) {
        for( pp in 2:p ) {
            glmformula <- paste( glmformula, " + I(x^", pp,")", sep = "");
        }
    }
    fit <- NULL;

# GLM fit
    if( link != "logit"      &&
            link != "probit"     &&
            link != "loglog"     &&
            link != "comploglog" &&
            link != "weibull"    &&
            link != "revweibull" ) {
            	linkfun <- link
            	}else{
            		linkfun <- paste( link, "_link_private", sep = "" );
            		}

    suppressWarnings( lapsing <- optim( initval, likfun )$par );
    lapsing <- 1 / ( 1 + exp( -lapsing ) );
    guessing <- 0;
    if( linkfun != "weibull_link_private" && linkfun != "revweibull_link_private" ) {
        tmpglm <- glm( glmformula, data = glmdata, weights = m,
                 family = binomial( eval( call( linkfun, guessing, lapsing ) ) ) );
    }
    else {
        tmpglm <- glm( glmformula, data = glmdata, weights = m,
                 family = binomial( eval( call( linkfun, K, guessing, lapsing ) ) ) );
    }
    
	b <- tmpglm$coeff
    
    tmpglm$df.residual <- length(x) - (p + 1) - 1
    
    fit <- tmpglm
    
    value <- NULL
    value$lapsing <- lapsing
    value$b <- b
    value$fit <- fit
    
    return( value );
}