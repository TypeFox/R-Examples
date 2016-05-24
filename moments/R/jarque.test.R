jarque.test <- function(x)
{
    if ( !is.vector( x ) )
        stop( "argument x is not a vector" )
    if ( !is.numeric( x ) )
        stop( "argument x is not numeric" )
    DNAME <- deparse( substitute( x ) )
    n <- length(x)
    ALTERNATIVE <- "greater"
    METHOD <- "Jarque-Bera Normality Test"
    K <- kurtosis( x )
    S <- skewness( x )
    JB  <- ( n / 6 ) * ( S^2 + 0.25 * ( ( K - 3 )^2 ) )
    pval <- 1 - pchisq( JB, df=2 )
    JBVAL <- list( statistic=c(JB=JB), p.value=pval, alternative=ALTERNATIVE, method=METHOD,
                   data.name=DNAME )
    class( JBVAL ) <- "htest"
    return( JBVAL )
}
   