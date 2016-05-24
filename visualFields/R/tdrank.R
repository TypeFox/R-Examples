tdrank <- function( td ) {
# gets the total-deviation rank curve

  tdr <- td
  for( i in 1:nrow( td ) ) {
    tdr[i,visualFields::vfsettings$locini:ncol( td )] <- td[i,( order( td[i,visualFields::vfsettings$locini:ncol( td )], decreasing = TRUE ) + visualFields::vfsettings$locini - 1 )]
  }
# Remove all columns that are NAs because they correspond to locations of the blind spot.
# This is a bit clumsy, but necessary if we allow vf-objects to have different patterns.
# That is, if some rows are for 24d2 and others for 30d2, for which the number of tested
# locations are different...
  for( i in ncol( tdr ):visualFields::vfsettings$locini ) {
    if( all( is.na( tdr[,i] ) ) ){
      tdr[,i] <- NULL
    } else {
      break
    }
  }

  return( tdr )
}
