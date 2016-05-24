tdpmap <- function( td ) {
# gets the probability score from the current normative value reference.

# init
  tdp <- td

# NOTE: if the cutoffs defined in vfenv$nv are e.g. 0.5, 1.0, 5.0, 95.0, and 100.0,
# a cutoff value of 0.5% means that the p-value at that location was smaller
# than 0.5%, a value of 1.0% means that the p-value of the location was between
# 0.5 and 1.0%, etc. A value of 95.0% means within normal limits and a value of
# 100.0 means above 95th percentile, i.e. location was "too sensitive"
  
  for( i in 1:nrow( tdp ) ) {
# get how many locations we need to look at
    texteval <- paste( "vfsettings$", tdp$tpattern[i], "$locnum", sep = "" )
    locnum <- eval( parse( text = texteval ) )
# get the reference values for the TD map...
    texteval <- paste( "vfenv$nv$", tdp$tpattern[i], "_", tdp$talgorithm[i], "$TDpercloc", sep = "" )
    tdco <- eval( parse( text = texteval ) )
# init TD probability maps
    tdp[i,visualFields::vfsettings$locini:( visualFields::vfsettings$locini - 1 + locnum )] <- NA
# lets start comparing values
    td_iter <- as.numeric( td[i,visualFields::vfsettings$locini:( visualFields::vfsettings$locini - 1 + locnum )] )
    idx <- which( td_iter <= tdco[,1] ) + visualFields::vfsettings$locini - 1
    if( length( idx ) > 0 ) {
      tdp[i,idx] <- visualFields::vfenv$nv$pmapsettings$cutoffs[1]
    }
    for( j in 2:( length( visualFields::vfenv$nv$pmapsettings$cutoffs ) - 1 ) ) {
      idx <- which( tdco[,j-1] < td_iter & td_iter <= tdco[,j] ) + visualFields::vfsettings$locini - 1
      if( length( idx ) > 0 ) {
        tdp[i,idx] <- visualFields::vfenv$nv$pmapsettings$cutoffs[j]
      }
    }
    idx <- which( td_iter > tdco[,j] ) + visualFields::vfsettings$locini - 1
    if( length( idx ) > 0 ) {
      tdp[i,idx] <- visualFields::vfenv$nv$pmapsettings$cutoffs[length( visualFields::vfenv$nv$pmapsettings$cutoffs )]
    }
  }

  return( tdp )
}
