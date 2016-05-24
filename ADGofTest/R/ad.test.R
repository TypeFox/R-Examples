#################################################################
# 
# File:         ad.test.r
# Purpose:      Implements the Anderson Darling GoF test
#
# Created:      20090625
# Author:       Carlos J. Gil Bellosta
#
# Modifications: 
#
#################################################################

ad.test <- function( x, distr.fun, ... )
{
    DNAME  <- deparse( substitute( x ) )

    x <- sort( x )

    if( missing( distr.fun ) && ( x[1] < 0 || x[ length( x ) ] > 1 ) )
        stop( paste( "Data ", DNAME, " is not in the [0,1] range."  ) )

    if( ! missing( distr.fun ) ){
        x <- distr.fun( x, ... )
        DNAME <- paste( DNAME, " and ", deparse( substitute( distr.fun ) ) ) 
    }

    STATISTIC <- ad.test.statistic( x )
    names( STATISTIC ) <- "AD"

    PVAL      <- 1 - ad.test.pvalue( STATISTIC, length( x ) )

    METHOD <- "Anderson-Darling GoF Test"
    ALTERNATIVE <- "NA"

    RVAL <- list(statistic = STATISTIC, p.value = PVAL, alternative = ALTERNATIVE, method = METHOD, data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)

}
