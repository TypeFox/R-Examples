
cv.svd.gabriel <- function(x, krow=2, kcol=2, 
                           maxrank = floor(min(n - n/krow, p - p/kcol))) {
    x <- as.matrix( x )
    n <- nrow(x)
    p <- ncol(x)
    
    if (n < 2)
        stop ("x should have at least two rows")
    if (p < 2)
        stop ("x should have at least two columns")
    if ((krow > n) || (krow <= 1)) 
        stop("krow outside allowable range")
    if ((kcol > p) || (kcol <= 1))
        stop("kcol outside allowable range")
    if (maxrank < 0)
        stop("maxrank should be non-negative")
    
    krow.o <- krow; krow <- round.fold(n, krow);
    kcol.o <- kcol; kcol <- round.fold(p, kcol);
    
    if (krow != krow.o) 
        warning("krow has been set to ", krow)
    if (kcol != kcol.o)
        warning("kcol has been set to ", kcol)
    
    s.r <- choose.sets(n, krow)
    s.c <- choose.sets(p, kcol)

    n0 <- n - max( table( s.r ) )
    p0 <- p - max( table( s.c ) )
    maxrank.o <- maxrank
    maxrank   <- min( n0, p0, round( maxrank.o ) )
    
    if (!missing(maxrank) && maxrank != maxrank.o)
        warning("maxrank has been set to ", maxrank)
    
    msept <- .Call( "R_cv_svd_gabriel", x, krow, kcol, maxrank, 
                     as.integer( s.r-1 ), as.integer( s.c-1 ) )
    msep  <- t(msept)
    colnames( msep ) <- 0:maxrank
    
    res    <- list( call=match.call(), krow=krow, kcol=kcol, maxrank=maxrank, 
                    msep=msep, rowsets=s.r, colsets=s.c)
    class( res ) <- c("cvsvd_gabriel", "cvsvd")
    res
}
