tdrankglm <- function( tdr, familytxt = c( "gaussian" ), link = make.link( "logit" ), rankCentral = NULL, scaleFactor = 52.4 ) {
# gets a generalized linear model fit to total-deviation curve

  tdrglm <- NULL

  tdrtofit <- tdr
  if( is.null( rankCentral ) ) {
    rankCentral <- c( 1:length( tdrtofit ) )
  } else{
    tdrtofit <- tdrtofit[rankCentral]
  }
  rankval <- c( 1:length( tdr ) )
# prepare to get glm fit
  glmfit   <- NULL
  glmfit$x <- -tdrtofit
  glmfit$y <- rankCentral / scaleFactor
  glmfit   <- as.data.frame( glmfit )
# formula
  formula <- y ~ x
# get the family for the 
  texteval       <- paste( familytxt, "( link = link )", sep = "" )
  family         <- eval( parse( text = texteval ) )
  glmodel        <- glm( formula, family = family, data = glmfit )
  tdrglm$glmodel <- glmodel
# get the fitted values
  tdrglm$rank <- rankval
  tdrglm$shape  <- as.numeric( - 1 / glmodel$coefficient[2] )
  tdrglm$height <- as.numeric( glmodel$coefficient[1] / glmodel$coefficient[2] )
  tdrglm$scale  <- scaleFactor
#  tdrglm$val  <- - ( link$linkfun( rankval / scaleFactor ) - glmodel$coefficients[1] ) / glmodel$coefficients[2]
  tdrglm$val    <- tdrglm$shape * link$linkfun( rankval / scaleFactor ) + tdrglm$height
# get rest of details
  tdrglm$link   <- link
  return( tdrglm )
}
