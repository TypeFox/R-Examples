ghpostd <- function( td, correction = FALSE, rankRef = c( "default" ) ) {

  if( nrow( td ) > 1) {
    stop( "pass only one visual field here" )
  }
# get blind-spot position
  texteval <- paste( "vfsettings$", td$tpattern, "$bs", sep = "" )
  bspos <- eval( parse( text = texteval ) )
# get how many locations we need to look at
  texteval <- paste( "vfsettings$", td$tpattern, "$locnum", sep = "" )
  locnum <- eval( parse( text = texteval ) ) - length( bspos )
  if( rankRef == "default" ) {
# find, for the pattern used which is the rank position corresponding,
# approximately (although not always quite) with the 85th TD percentile
    texteval <- paste( "vfsettings$", td$tpattern, "$locrPD", sep = "" )
    rankRef <- eval( parse( text = texteval ) )
    if( is.null( rankRef ) ) {
      stop( "pattern of stimuli locations not recognized" )
    }
  }
# get TD values in a list and remove blind spot
  tdr <- tdrank( td )
# rankRef-th largest TD value
  gh <- as.numeric( tdr[visualFields::vfsettings$locini + rankRef - 1] )
  if( correction ) {
    texteval <- paste( "vfenv$nv$", td$tpattern, "_", td$talgorithm, "$nvtdrank$mtdr", sep = "" )
    tdrref  <- eval( parse( text = texteval ) )
    gh <- tdrref[rankRef] - gh
  }
  return( gh )
}
