pdpmapghr <- function( pd ) {
# gets the probability score from the current normative value reference.

# init
  pdp <- pd

# NOTE: if the cutoffs defined in vfenv$nv are e.g. 0.5, 1.0, 5.0, 95.0, and 100.0,
# a cutoff value of 0.5% means that the p-value at that location was smaller
# than 0.5%, a value of 1.0% means that the p-value of the location was between
# 0.5 and 1.0%, etc. A value of 95.0% means within normal limits and a value of
# 100.0 means above 95th percentile, i.e. location was "too sensitive"
  
  for( i in 1:nrow( pdp ) ) {
# get how many locations we need to look at
    texteval <- paste( "vfsettings$", pdp$tpattern[i], "$locnum", sep = "" )
    locnum <- eval( parse( text = texteval ) )
# get the reference values for the PD map...
    texteval <- paste( "vfenv$nv$", pdp$tpattern[i], "_", pdp$talgorithm[i], "$PDGHRpercloc", sep = "" )
    pdco <- eval( parse( text = texteval ) )
# init PD probability maps
    pdp[i,visualFields::vfsettings$locini:( visualFields::vfsettings$locini - 1 + locnum )] <- NA
# lets start comparing values
    pd_iter <- as.numeric( pd[i,visualFields::vfsettings$locini:( visualFields::vfsettings$locini - 1 + locnum )] )
    idx <- which( pd_iter <= pdco[,1] ) + visualFields::vfsettings$locini - 1
    if( length( idx ) > 0 ) {
      pdp[i,idx] <- visualFields::vfenv$nv$pmapsettings$cutoffs[1]
    }
    for( j in 2:( length( visualFields::vfenv$nv$pmapsettings$cutoffs ) - 1 ) ) {
      idx <- which( pdco[,j-1] < pd_iter & pd_iter <= pdco[,j] ) + visualFields::vfsettings$locini - 1
      if( length( idx ) > 0 ) {
        pdp[i,idx] <- visualFields::vfenv$nv$pmapsettings$cutoffs[j]
      }
    }
    idx <- which( pd_iter > pdco[,j] ) + visualFields::vfsettings$locini - 1
    if( length( idx ) > 0 ) {
      pdp[i,idx] <- visualFields::vfenv$nv$pmapsettings$cutoffs[length( visualFields::vfenv$nv$pmapsettings$cutoffs )]
    }
  }

  return( pdp )
}