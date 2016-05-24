obfuscate <- function( obj=NULL ) {
  if( is.phe( obj ) ) {
    ## can't be a symbolic object here
    if( is.sym( obj ) )
      obj <- read.phe( get.sym(obj), sym=FALSE );
    
    ## obfuscate the names
    llen <- ncol(obj);
    names(obj)[3:llen] <- as.character(floor(runif(llen-2)*10000));

    ## obfuscate the data
    for( col in 3:llen ) {
      obj[,col] <- sample( obj[,col] );
    }

    return( obj );
  }else if( is.ped( obj ) || is.pedlist(obj) ){
    if( is.pped( obj ) )
      stop( "pped objects not yet supported, only ped objects (the decompressed version." );

    ## can't be a symbolic object since we're modifying it!
    if( is.sym( obj ) )
      obj <- read.ped( get.sym(obj), sym=FALSE );
    if( is.pedlist( ped ) )
      ped <- as.ped( ped );
    
    ## obfuscate the names
    llen <- ncol(obj);
    numsnp <- (llen-6)/2;
    marker.names <- paste( "m", as.character(floor(runif((llen-6)/2)*10000)), sep="" );
    col.names <- names(obj)[1:6]  ## the immutable names
    for( i in 1:numsnp ){
      col.names <- c(col.names,
                     paste(marker.names[i], ".a", sep=""),
                     paste(marker.names[i], ".b", sep=""));
    }
    names( obj ) <- col.names;

    ## obfuscate the data
    for( col in 7:llen ) {
      obj[,col] <- sample( obj[,col] );
    }

    return( obj );
  }

  stop( "Run with either a 'phe' or 'ped' object, but ensure that they are loaded with the option 'sym=FALSE'." );
}
