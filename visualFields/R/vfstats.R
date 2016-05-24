vfstats <- function( vf ) {
# obtains statistical global indices for visual fields. The returned object has the same info columns
# as vf + msens, ssens, mtdev, stdev, mpdev, spdev for mean and standard deviations for sensitivities,
# TD values and PD values

  numstatindices <- 6 # change if more are to be included here

# init objects
  vfs    <- vf[,1:( visualFields::vfsettings$locini-1 )]
  vfsaux <- NULL

# get TD and PD values
  td <- tdval( vf )
  pd <- pdval( td )

# for each row, we need to find the corresponding standard deviations for the 
  for( i in 1:nrow( vfs ) ) {
# get how many locations we need to look at
    texteval <- paste( "vfsettings$", td$tpattern[i], "$locnum", sep = "" )
    locnum <- eval( parse( text = texteval ) )
# get blind spot
    texteval <- paste( "vfsettings$", td$tpattern[i], "$bs", sep = "" )
    bs <- eval( parse( text = texteval ) )
# get weights
    texteval <- paste( "vfenv$nv$", td$tpattern[i], "_", td$talgorithm[i], "$sds", sep = "" )
    wgt <- 1 / eval( parse( text = texteval ) )
# senect sensitivity, TD and PD values for the visual field in this iteration
    vf_iter <- as.numeric(vf[i,visualFields::vfsettings$locini:( visualFields::vfsettings$locini - 1 + locnum )] )
    td_iter <- as.numeric(td[i,visualFields::vfsettings$locini:( visualFields::vfsettings$locini - 1 + locnum )] )
    pd_iter <- as.numeric(pd[i,visualFields::vfsettings$locini:( visualFields::vfsettings$locini - 1 + locnum )] )
# remove BS from everywhere
    if( all( !is.na( bs ) ) ) {
      wgt     <- wgt[-bs,]
      vf_iter <- vf_iter[-bs]
      td_iter <- td_iter[-bs]
      pd_iter <- pd_iter[-bs]
    }
# finally, get mean and std of sensitivity...
    vfsaux$msens[i] <- weighted.mean( vf_iter, w = wgt$sens )
    vfsaux$ssens[i] <- sqrt( wtd.var( vf_iter, weights = wgt$sens, normwt = TRUE ) )
# ... mean and std of TD values ...
    vfsaux$mtdev[i] <- weighted.mean( td_iter, w = wgt$td )
    vfsaux$stdev[i] <- sqrt( wtd.var( td_iter, weights = wgt$td, normwt = TRUE ) )
# ... and mean and std of PD values
    vfsaux$mpdev[i] <- weighted.mean( pd_iter, w = wgt$pd )
    vfsaux$spdev[i] <- sqrt( wtd.var( pd_iter, weights = wgt$pd, normwt = TRUE ) )
  }
  
  vfsaux <- as.data.frame( vfsaux )
  vfs[,visualFields::vfsettings$locini:( visualFields::vfsettings$locini - 1 + numstatindices )] <- vfsaux

  return( vfs )

}