"panel.dice" <-
function(x,y){
  tmp.cols <- c("Red","Green","Blue","Black","Yellow",
                "Purple","Orange","Brown","Grey","White")
  box.x <- c( 0.1, 0.9, 0.9, 0.1, 0.1 )
  box.y <- c( 0.1, 0.1, 0.9, 0.9, 0.1 )
  pips.x <- c( 0.5, 0.3, 0.7, 0.3, 0.7, 0.3, 0.7 )
  pips.y <- c( 0.5, 0.7, 0.3, 0.3, 0.7, 0.5, 0.5 )
  xx <- ceiling(sqrt(length(x)))
  yy <- ceiling( length(x)/xx )
  for( i in seq(along=x) ){
    xo <- y[i] %% xx
    yo <- yy-1-(y[i] %/% xx)
    lattice::llines( box.x+xo, box.y+yo,col=tmp.cols[i] )
    which <- c( x[i] %%2 == 1, x[i] > c(1,1,3,3,5,5) )
    tmp.x <- pips.x[which]
    tmp.y <- pips.y[which]
    if( runif(1) < 0.5 ) {
      tmp.x <- 1-tmp.x
    }
    if( runif(1) < 0.5 ) {
      tmp <- tmp.x
      tmp.x <- tmp.y
      tmp.y <- tmp
    }
    lattice::lpoints( tmp.x+xo, tmp.y+yo, pch=16,col='black')
  }
}

