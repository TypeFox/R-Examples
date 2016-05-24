vfstatspmap <- function( vfindices ) {
# gets the probability maps for each global index

  vfindicesp <- vfindices

# for each row, we need to find the corresponding standard deviations for the 
  for( i in 1:nrow( vfindicesp ) ) {
# get the reference values for the global-indices map...
    texteval <- paste( "vfenv$nv$", vfindicesp$tpattern[i], "_", vfindicesp$talgorithm[i], "$percglo", sep = "" )
    gico <- eval( parse( text = texteval ) )
    vfindicesp_iter <- as.numeric( vfindicesp[i,visualFields::vfsettings$locini:( visualFields::vfsettings$locini - 1 + nrow( gico ) )] )

# analysis of mean and std is different. This patch takes that into account.
# this is very clumsy way of doing it and could be improved
    mcols = which( strtrim( names( vfindicesp[i,visualFields::vfsettings$locini:( visualFields::vfsettings$locini - 1 + nrow( gico ) )] ), 1 ) == "m" )
    scols = which( strtrim( names( vfindicesp[i,visualFields::vfsettings$locini:( visualFields::vfsettings$locini - 1 + nrow( gico ) )] ), 1 ) == "s" )
    idx <- mcols[which( vfindicesp_iter[mcols] <= gico[mcols,1] )] + visualFields::vfsettings$locini - 1
    if( length( idx ) > 0 ) {
      vfindicesp[i,idx] <- visualFields::vfenv$nv$globalco[1]
    }
    idx <- scols[which( vfindicesp_iter[scols] >= gico[scols,1] )] + visualFields::vfsettings$locini - 1
    if( length( idx ) > 0 ) {
      vfindicesp[i,idx] <- visualFields::vfenv$nv$globalco[1]
    }
    for( j in 2:( length( visualFields::vfenv$nv$globalco ) - 1 ) ) {
      idx <- mcols[which( gico[mcols,j-1] < vfindicesp_iter[mcols] & vfindicesp_iter[mcols] <= gico[mcols,j] )] + visualFields::vfsettings$locini - 1
      if( length( idx ) > 0 ) {
        vfindicesp[i,idx] <- visualFields::vfenv$nv$globalco[j]
      }
      idx <- scols[which( gico[scols,j-1] >= vfindicesp_iter[scols] & vfindicesp_iter[scols] > gico[scols,j] )] + visualFields::vfsettings$locini - 1
      if( length( idx ) > 0 ) {
        vfindicesp[i,idx] <- visualFields::vfenv$nv$globalco[j]
      }
    }
    idx <- mcols[which( vfindicesp_iter[mcols] > gico[mcols,j] )] + visualFields::vfsettings$locini - 1
    if( length( idx ) > 0 ) {
      vfindicesp[i,idx] <- visualFields::vfenv$nv$globalco[length( visualFields::vfenv$nv$globalco )]
    }
    idx <- scols[which( vfindicesp_iter[scols] < gico[scols,j] )] + visualFields::vfsettings$locini - 1
    if( length( idx ) > 0 ) {
      vfindicesp[i,idx] <- visualFields::vfenv$nv$globalco[length( visualFields::vfenv$nv$globalco )]
    }
  }

  return( vfindicesp )

}
