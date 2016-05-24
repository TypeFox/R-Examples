#$Id: genomedistances.R 176 2014-08-04 11:15:52Z larssn $



fluidity <- function( pan.matrix, n.sim=10 ){
  pan.matrix[which( pan.matrix > 0, arr.ind=T )] <- 1
  ng <- dim( pan.matrix )[1]
  flu <- rep( 0, n.sim )
  for( i in 1:n.sim ){
    ii <- sample( ng, 2 )
    g1 <- pan.matrix[ii[1],]
    g2 <- pan.matrix[ii[2],]
    flu[i] <- (sum( g1>0 & g2==0 )+sum( g1==0 & g2>0 ))/(sum(g1)+sum(g2))
  }
  flu.list <- list( Mean=mean( flu ), Std=sd( flu ) )
  return( flu.list )
}


distJaccard <- function( pan.matrix ){
  pan.matrix[which( pan.matrix > 0, arr.ind=T )] <- 1
  ng <- dim( pan.matrix )[1]
  dtab <- matrix( 0, nrow=ng, ncol=ng )
  rownames( dtab ) <- rownames( pan.matrix )
  colnames( dtab ) <- rownames( pan.matrix )
  for( i in 1:(ng-1)){
    g1 <- pan.matrix[i,]
    for( j in (i+1):ng ){
      cs <- g1+pan.matrix[j,]
      dtab[j,i] <- dtab[i,j] <- 1 - sum( cs>1 )/sum( cs>0 )
    }
  }
  return( as.dist( dtab ) )
}


distManhattan <- function( pan.matrix, scale=0.0, weights=rep( 1, dim( pan.matrix )[2] ) ){
  if( (scale>1)|(scale<0) ){
    warning( "scale should be between 0.0 and 1.0, using scale=0.0" )
    scale <- 0.0
  }
  idx <- which( pan.matrix > 0, arr.ind=T )
  pan.matrix[idx] <- 1 + (pan.matrix[idx]-1)*scale
  
  pan.matrix <- pan.matrix * matrix( weights, nrow=dim( pan.matrix )[1], ncol=dim( pan.matrix )[2], byrow=T )
  dtab <- dist( pan.matrix, method="manhattan" )
  return( dtab )
}


geneWeights <- function( pan.matrix, type=c("shell","cloud") ){
  ng <- dim( pan.matrix )[1]
  nf <- dim( pan.matrix )[2]
  pan.matrix[which( pan.matrix>0, arr.ind=T )] <- 1
  cs <- colSums( pan.matrix )
  
  midx <- grep( type[1], c( "shell", "cloud" ) )
  if( length( midx ) == 0 ){
    warning( "Unknown weighting:", type, ", using shell weights" )
    midx <- 1
  }
  W <- rep( 1, nf )
  x <- 1:ng
  ww <- 1/(1+exp( ((x-1)-(max(x)-1)/2)/((max(x)-1)/10) ))
  if( midx == 1 ) ww <- 1-ww
  for( i in 1:ng ) W[which( cs == i )] <- ww[i]
  return( W )
}
