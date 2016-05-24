lidLensArtifact <- function( vf, min_dB = c( 12 ) ) {
# I wish I knew how to do this well... and this approach needs to be improved as
# it is just too very coarse. ONLY TO BE USED WITH CONTROL DATA!!!

  idx_lidlens <- NULL
  k <- 0
  for( i in 1:nrow( vf ) ) {
# remove blind-spot data
    texteval <- paste( "vfsettings$", vf$tpattern[i], "$bs", sep = "" )
    bs <- eval( parse( text = texteval ) )
    vf2 <- vf[i,visualFields::vfsettings$locini:ncol( vf )]
    vf2 <- vf2[,-bs]
    if( length( which( vf2 <= min_dB ) ) >= 1 ) {
      k <- k + 1
      idx_lidlens[k] <- i
    }
  }

  return( idx_lidlens )
}
