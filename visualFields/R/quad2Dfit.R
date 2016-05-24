quad2Dfit <- function( val, patternMap = visualFields::saplocmap$p24d2,
                       bspos = visualFields::vfsettings$p24d2$bs ) {
# fit a 2D quadratic function using values in val as "observations" for the x and y coordinates in patternMap

  val           <- as.data.frame( val )
  val_sm        <- NULL
  val_sm$result <- val$val
  val_sm$x      <- patternMap$xod
  val_sm$y      <- patternMap$yod
  val_sm        <- as.data.frame( val_sm )

# remove blind spot
  if( all( !is.na( bspos ) ) ) val_sm <- val_sm[-bspos,]

  formula <- result ~ x + y + I( x^2 ) + I( y^2 )
#  val_sm$result <- predict( lm( formula, data = val_sm ) )

# update new values for slopes and intercepts witht the "smoothed" ones
  val[attr( val_sm, "row.names" ),] <- predict( lm( formula, data = val_sm ) )

  return( c( val$val ) )
}
