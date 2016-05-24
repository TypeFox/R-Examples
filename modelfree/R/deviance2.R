deviance2 <- function( r, m, pfit ) {
#
# The function calculates the deviance for the fitted values of the psychometric function pfit.
#  
#
# INPUT
# r    - number of successes
# m    - number of trials
# pfit - fittd values
#
# OUTPUT
# D - deviance

# Both arguments are mandatory
    if( missing("pfit") || missing("r") || missing("m") ) {
        stop("Check input. First 3 arguments are mandatory");
    }

# adjustment to avoid degenerate values
    r[which(r >= m)]<-r[which(r >= m)] - .001;
    r[which(r <= 0)]<-.001;

    pfit[which(pfit >= 1)]<-1 - .001;
    pfit[which(pfit <= 0)]<-.001;   

# deviance
    return(2 * sum( ( r * log( r / ( m * pfit) ) + ( m - r ) * log( ( m - r ) / ( m - m * pfit ) ) ) ));
}
