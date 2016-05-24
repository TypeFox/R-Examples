vfsegmentcoord <- function( lineMap, length = 2.5 ) {

  length <-  length / 2

  coords <- NULL
  coords$x1[1:nrow( lineMap )] <- c( NA )
  coords$y1[1:nrow( lineMap )] <- c( NA )
  coords$x2[1:nrow( lineMap )] <- c( NA )
  coords$y2[1:nrow( lineMap )] <- c( NA )
  coords                       <- as.data.frame( coords )

  coords$x1 <- length / sqrt( 1 + ( lineMap$jmslope )^2 ) + lineMap$xod
  coords$y1 <- lineMap$jmslope * ( length / sqrt( 1 + ( lineMap$jmslope )^2 ) ) + lineMap$yod
  coords$x2 <- 2 * lineMap$xod - coords$x1
  coords$y2 <- 2 * lineMap$yod - coords$y1
  
  return ( coords )
}