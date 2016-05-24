vfindex <- function( vf, td2pdcutoff = -20, perc = 5, vfiset = visualFields::vfidefault ) {
# calculates the vfi. It is flexible enought to change the weithgs, weighting regions,
# the cutoff value for MD and, even, the percentile to look at whether the 85th TD
# percentile is within normal limits or not. It returns the mean AND the SD of vfi
# obtained for each location

  numstatindices <- 2 # change if more are to be included here

# init
  vfi    <- vf[,1:( visualFields::vfsettings$locini-1 )]
  vfiaux <- NULL
# calculate TD and PD maps
  td <- tdval( vf )
  pd <- pdval( td )
# get TD and PD probability maps
  tdp <- tdpmap( td )
  pdp <- pdpmap( pd )
# calculate global vf indices
  vfindices <- vfstats( vf )
  for( i in 1:nrow( vf ) ) {
# get how many locations we need to look at
    texteval <- paste( "vfsettings$", vf$tpattern[i], "$locnum", sep = "" )
    locnum <- eval( parse( text = texteval ) )
# init
    vfiloc <- numeric( locnum )
    wgt    <- numeric( locnum )
# get blind spot
    texteval <- paste( "vfsettings$", vf$tpattern[i], "$bs", sep = "" )
    bs <- eval( parse( text = texteval ) )
# find, for the pattern used which is the rank position corresponding,
# approximately (although not always quite) with the 85th TD percentile
    texteval <- paste( "vfsettings$", pd$tpattern[i], "$locrPD", sep = "" )
    rankRef <- eval( parse( text = texteval ) )
# get the norm data and calculate normal age-corrected sensitivities
    texteval <- paste( "vfenv$nv$", td$tpattern[i], "_", td$talgorithm[i], "$agelm", sep = "" )
    agelm <- eval( parse( text = texteval ) )
# get weights
    texteval <- paste( "vfiset$", td$tpattern[i], sep = "" )
    wtdaux <- eval( parse( text = texteval ) )
    for( j in 1:nrow( wtdaux$regweights ) ) {
      wgt[which( wtdaux$locregions$region == wtdaux$regweights$region[j] )] <- wtdaux$regweights$weight[j]
    }
# get sensitivities, td and pd values, and td and pd probability maps for the vf
# of this iteration
    vf_iter  <- as.numeric(vf[i,visualFields::vfsettings$locini:( visualFields::vfsettings$locini - 1 + locnum )] )
    td_iter  <- as.numeric(td[i,visualFields::vfsettings$locini:( visualFields::vfsettings$locini - 1 + locnum )] )
    pd_iter  <- as.numeric(pd[i,visualFields::vfsettings$locini:( visualFields::vfsettings$locini - 1 + locnum )] )
    tdp_iter <- as.numeric(tdp[i,visualFields::vfsettings$locini:( visualFields::vfsettings$locini - 1 + locnum )] )
    pdp_iter <- as.numeric(pdp[i,visualFields::vfsettings$locini:( visualFields::vfsettings$locini - 1 + locnum )] )
# remove blind spot from everywhere
    if( all( !is.na( bs ) ) ) {
      vf_iter  <- vf_iter[-bs]
      td_iter  <- td_iter[-bs]
      pd_iter  <- pd_iter[-bs]
      tdp_iter <- tdp_iter[-bs]
      pdp_iter <- pdp_iter[-bs]
      agelm    <- agelm[-bs,]
      wgt      <- wgt[-bs]
      vfiloc   <- vfiloc[-bs]
    }
# get age-corrected norms
    vf_age <- agelm$intercept + agelm$slope * vf$sage[i]
# check if we should use TD or PD values
    if( !( vfindices$mtdev[i] < td2pdcutoff & tdp_iter[order( td_iter, decreasing = TRUE )[rankRef]] <= perc ) ) {
# use PD values
      idx <- which( pdp_iter > perc )
      vfiloc[idx] <- 100
      idx <- which( pdp_iter <= perc )
      vfiloc[idx] <- 100 * ( 1 - abs( td_iter[idx] ) / vf_age[idx] )
      idx <- which( vf_iter == 0 )
      vfiloc[idx] <- 0
      } else {
# use TD values
      idx <- which( tdp_iter > perc )
      vfiloc[idx] <- 100
      idx <- which( tdp_iter <= perc )
      vfiloc[idx] <- 100 * ( 1 - abs( td_iter[idx] ) / vf_age[idx] )
      idx <- which( vf_iter == 0 )
      vfiloc[idx] <- 0
    }
    vfiaux$mvfi[i] <- weighted.mean( vfiloc, w = wgt )
    vfiaux$svfi[i] <- sqrt( wtd.var( vfiloc, weights = wgt, normwt = TRUE ) )
  }

  vfiaux <- as.data.frame( vfiaux )
  vfi[,visualFields::vfsettings$locini:( visualFields::vfsettings$locini - 1 + numstatindices )] <- vfiaux

  return( vfi )
}
