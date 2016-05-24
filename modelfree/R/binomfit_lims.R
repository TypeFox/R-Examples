binomfit_lims <- function( r, m, x, p = 1, link = c( "logit" ), guessing = 0, lapsing = 0,
                           K = 2 ) {
#
# The function fits a binomial generalised liner model with fixed guessing and lapsing rates.
#
# INPUT
#
# r    - number of successes at points x
# m    - number of trials at points x 
# x    - stimulus levels
#
# OPTIONAL INPUT
#
# p    - degree of the polynomial; default is p = 1 
# link - name of the link function; default is "logit"
# guessing - guessing rate; default is 0
# lapsing  - lapsing rate; default is 0
# K    - power parameter for Weibull and reverse Weibull link; default is 2
#
# OUTPUT
# 
# Object with 2 components: 
# b - vector of estiamted coefficients for the linear part 
# fit - glm object to be used in evaluation of fitted values

# MAIN PROGRAM
# First 3 arguments are mandatory
    # First 3 arguments are mandatory
    if( missing("x") || missing("r") || missing("m") ) {
        stop("Check input. First 3 arguments are mandatory");
    }

# CHECK ROBUSTNESS OF INPUT PARAMETERS
    checkdata<-list();
    checkdata[[1]] <- x;
    checkdata[[2]] <- r
    checkdata[[3]] <- m
    checkinput( "psychometricdata", checkdata );
    rm( checkdata )
    checkinput( "linkfunction", link );
    pn <- list()
	pn[[1]] <- p
	pn[[2]] <- x
    checkinput( "degreepolynomial", pn );
    checkinput( "guessingandlapsing", c( guessing, lapsing ) );
    if( link == "weibull"  || link == "revweibull") {
    	checkinput( 'exponentk', K );
    	}


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

    if( linkfun != "weibull_link_private" && linkfun != "revweibull_link_private" ) {
        fit <- glm( glmformula, data = glmdata, weights = m,
                    family = binomial( eval( call( linkfun, guessing, lapsing ) ) ) );
    }
    else {
        fit <- glm( glmformula, data = glmdata, weights = m,
                    family = binomial( eval( call( linkfun, K, guessing, lapsing ) ) ) );
    }
    
    value <- NULL       
    value$b <- fit$coeff
    value$fit <- fit
    
    return( value );
}
