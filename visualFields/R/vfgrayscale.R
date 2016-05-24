vfgrayscale <- function( sens, age, pattern, algorithm ) {

# get the norm data and calculate normal age-corrected sensitivities
  texteval <- paste( "vfenv$nv$", pattern, "_", algorithm, "$agelm", sep = "" )
  agelm    <- eval( parse( text = texteval ) )
  vals     <- agelm$intercept + agelm$slope * age
  idx      <- which( is.na( vals ) )
  if( length( idx ) > 0 ) vals <- vals[-idx]
  sens <- as.numeric( sens )
  sens <- sens / max( vals )
  sens[which( sens < 0 )] <- 0
  sens[which( sens > 1 )] <- 1
  return ( matrix( rep( sens, 3 ), nrow = length( sens ), ncol = 3 ) )

}
