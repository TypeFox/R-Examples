
meas.est <- function( datameas, id, data=NULL )
{
    if ( nargs() != 2 )
    {
        stop("An id vector is required, to identify which subject each measurement belongs to.")
    }
#     if ( is.null(data)==FALSE )
#     {
#         attach(data)
#     }
    
    
    if(!is.null(data))
      stop("'data' argument no longer supported.")
    
    
    datameas <- as.matrix( datameas )
    siz <- dim( datameas )
    if ( length(id)!=siz[1] )
    {
        stop("The id vector must have the same number of rows as the data matrix")
    }

    idlabels <- sort( unique( id ) )
    n        <- length( idlabels )
    ni       <- matrix( 0, n, siz[2] )
    dat      <- matrix( NA, n, siz[2] )
    vrs      <- rep( NA, siz[2] * siz[2] * n )
    dim(vrs) <- c( siz[2], siz[2], n )

    is.OK    <- is.finite( apply( datameas , 1, sum ) ) #the rows with no nan's

    for ( i in 1:n )
    {
        ref      <- id==idlabels[i] & is.OK
        ni[i]    <- sum( as.numeric(ref) )
        dat[i,]  <- apply( as.matrix( datameas[ ref, ] ), 2, mean )
        if ( ni[i] > 1 )
            { vrs[, , i] <- var( datameas[ ref, ] ) / ni[i] }
    }

    V <- apply(vrs, 1:2, mean, na.rm=TRUE)

#     if ( is.null(data)==FALSE )
#     {
#        detach(data)
#     }

    list( V=V, dat.mean=dat )
}
