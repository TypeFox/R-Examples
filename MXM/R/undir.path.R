############################
############################
##### Returns all nodes within all undirected path between any two nodes
#####
############################
############################

undir.path <- function(G, y, x) {

  ## G is an adjacency matrix
  ## y and x are two nodes

  geita <- geit <- which( G[x, ] == 1 )
  geita2 <- geit2 <- which( G[y, ] == 1 )
  
  if ( length(geita) == 0 & length(geita2) == 0 ) {
    geit <- c(y, x)
  } else if ( length(geita) == 0 & length(geita2) > 0 ) {
    geit <- rbind(c( y, geita2), c(x, numeric( length(geita2) ) ) )
  } else if ( length(geita) > 0 & length(geita2) == 0 ) {  
    geit <- rbind(c( x, geita), c(y, numeric( length(geita) ) ) )
  } else {

    for (i in geit) geita <- c( geita, which(G[i, ] == 1) )
    while ( length(intersect(y, geita) ) == 0 ) {
      for (i in geita) geita <- c( geita, which(G[i, ] == 1) )
      geita <- unique(geita)
    } 

    for (i in geit2) geita2 <- c( geita2, which(G[i, ] == 1) )
    while ( length(intersect(x, geita2) ) == 0 ) {
      for (i in geita2) geita2 <- c( geita2, which(G[i, ] == 1) )
      geita2 <- unique(geita2) 
    }

    geit <- intersect(geita, geita2)
    geit <- geit[ geit != y ]
    geit <- geit[ geit != x ]
    ina <- 1:length(geit)
 
    for ( i in geit ) {
      a <- which( G[i, ] == 1 )
      xin <- length( intersect(a, x) )
      yin <- length( intersect(a, y) )
      gein <- length( intersect(a, geit) )
      if ( xin + yin +  gein <= 1 ) {
        ina[ which(geit == i) ] <- 0
      }
    }
    geit <- c(y, geit[ ina ], x)
  }
  
  geit 
}