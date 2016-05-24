# wrapper function for outlier detection methods.
# 23.12.2009 version 1, mvdl
getOutliers <- function(y, method="I",  ...)
{
    if ( !(method %in% c("I","II") ) )
        stop("method not recognized (choose I or II)")
    out <- switch( method,
        I = getOutliersI(y, ...),
        II = getOutliersII(y, ...)
        )
    return(out)
}


