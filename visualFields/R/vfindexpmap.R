vfindexpmap <- function( vfi ) {
# gets the probability maps for each global index

  vfip <- vfi

  idxm <- visualFields::vfsettings$locini
  idxs <- visualFields::vfsettings$locini + 1

# for each row, we need to find the corresponding standard deviations for the 
  for( i in 1:nrow( vfip ) ) {
# get the reference values for the global-indices map...
    texteval <- paste( "vfenv$nv$", vfip$tpattern[i], "_", vfip$talgorithm[i], "$percvfi", sep = "" )
    gico <- eval( parse( text = texteval ) )
    vfip_iter <- as.numeric( vfip[i,visualFields::vfsettings$locini:( visualFields::vfsettings$locini - 1 + nrow( gico ) )] )

# analysis of mean and std is different. This patch takes that into account.
    if( vfip_iter[1] <= gico[1,1] ) {
      vfip[i,idxm] <- visualFields::vfenv$nv$globalco[1]
    }
    if( vfip_iter[2] >= gico[2,1] ) {
      vfip[i,idxs] <- visualFields::vfenv$nv$globalco[1]
    }
    for( j in 2:( length( visualFields::vfenv$nv$globalco ) - 1 ) ) {
      if( gico[1,j-1] < vfip_iter[1] & vfip_iter[1] <= gico[1,j] ) {
        vfip[i,idxm] <- visualFields::vfenv$nv$globalco[j]
      }
      if( gico[2,j-1] >= vfip_iter[2] & vfip_iter[2] > gico[2,j] ) {
        vfip[i,idxs] <- visualFields::vfenv$nv$globalco[j]
      }
    }
    if( vfip_iter[1] > gico[1,j] ) {
      vfip[i,idxm] <- visualFields::vfenv$nv$globalco[length( visualFields::vfenv$nv$globalco )]
    }
    if( vfip_iter[2] < gico[2,j] ) {
      vfip[i,idxs] <- visualFields::vfenv$nv$globalco[length( visualFields::vfenv$nv$globalco )]
    }
  }

  return( vfip )

}
