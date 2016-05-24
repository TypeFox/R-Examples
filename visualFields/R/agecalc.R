agecalc <- function( from, to, daysyear = NULL )  { 

  if( !is.null( daysyear ) ) {
    return( round( as.numeric( to - from ) / daysyear ) )
  }
  ltfrom <- as.POSIXlt( from )
  ltto   <- as.POSIXlt( to )
  
  age <- ltto$year - ltfrom$year
  # if month has not been reached yet, then the person is a year younger
  idx <- which( ltto$mon < ltfrom$mon )
  if( length( idx ) > 0 ) age[idx]  <- age[idx] - 1
  idx <- which( ltto$mon == ltfrom$mon & ltto$mday < ltfrom$mday )
  if( length( idx ) > 0 ) age[idx]  <- age[idx] - 1
  
  return( age )
}