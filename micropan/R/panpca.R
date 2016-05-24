#$Id: panpca.R 175 2014-07-27 15:57:47Z larssn $


panpca <- function( pan.matrix, scale=0.0, weights=rep( 1, dim( pan.matrix )[2] ) ){
  if( (scale>1)|(scale<0) ){
    warning( "scale should be between 0.0 and 1.0, using scale=0.0" )
    scale <- 0.0
  }
  idx <- which( pan.matrix > 0, arr.ind=T )
  pan.matrix[idx] <- 1 + (pan.matrix[idx]-1)*scale
  pan.matrix <- pan.matrix * matrix( weights, nrow=dim( pan.matrix )[1], ncol=dim( pan.matrix )[2], byrow=T )
  idx <- which( apply( pan.matrix, 2, sd ) > 0 )
  X <- pan.matrix[,idx]

  pca <- prcomp( X )
  evar <- pca$sdev^2/sum( pca$sdev^2 )
  scores <- pca$x
  loadings <- pca$rotation

  pan.pca <- list( Evar=evar, Scores=scores, Loadings=loadings, Scale=scale, Weights=weights )
  class( pan.pca ) <- c( "Panpca", "list" )
  return( pan.pca )
}

plot.Panpca <- function( x, cum=FALSE, col="black", ... ){
  Panpca <- x
  if( cum ){
    y <- cumsum( Panpca$Evar )
  } else {
    y <- Panpca$Evar
  }
  barplot( y, col=col, names.arg=1:length( y ), xlab="Principal components", ylab="Relative explained variance", ... )
}

summary.Panpca <- function( object, ... ){
  Panpca <- object
  for( i in 1:length( Panpca$Evar ) ){
    cat( "Principal component ", i, " explains ", round( 1000*Panpca$Evar[i] )/10, "% of the variation\n", sep="" )
  }
}

str.Panpca <- function( object, ... ){
  Panpca <- object
  cat( "Pan-matrix PCA over", dim( Panpca$Scores )[1], "genomes\n" )
}



plotScores <- function( pan.pca, x=1, y=2, show.labels=TRUE, labels=NULL, col="black", pch=16, ... ){
  Z <- pan.pca$Scores
  xr <- range( Z[,x] )
  yr <- range( Z[,y] )
  args <- list(...)
  ii <- match( c("xlab","ylab"), names( args ) )
  if( is.na(ii[1]) ){
    xlab=paste( "PC", x, " (", round( 100*pan.pca$Evar[x] ), "%)", sep="" )
  } else {
    xlab <- args[[ii[1]]]
  }
  if( is.na(ii[2]) ){
    ylab=paste( "PC", y, "(", round( 100*pan.pca$Evar[y] ), "%)", sep="" )
  } else {
    ylab <- args[[ii[2]]]
  }
  
  plot( xr, c(0,0), type="l", col="gray", xlim=xr, ylim=yr, xlab=xlab, ylab=ylab )
  points( c(0,0), yr, type="l", col="gray" )
  
  pm.gid <- rownames( Z )
  if( length( col )>1 ){
    if( is.null( names( col ) ) ) stop( "Each element in col must be named by its GID.tag" )
    gid <- names( col )
    idx <- match( pm.gid, gid )
    if( sum( is.na( idx ) ) > 0 ) stop( "GID.tag names does not match GID.tags in the pan-matrix" )
    cols <- col[idx]
    if( length( cols ) != length( pm.gid ) ) stop( "The number of elements in col does not correspond to the number of genomes in the pan-matrix" )
  } else {
    cols <- rep( col, length.out=length( pm.gid ) )
  }
  if( is.null( labels ) ){
    labs <- rownames( Z )
  } else {
    if( is.null( names( labels ) ) ) stop( "Each element in labels must be named by its GID.tag" )
    gid <- names( labels )
    idx <- match( pm.gid, gid )
    if( sum( is.na( idx ) ) > 0 ) stop( "GID.tag names does not match GID.tags in the pan-matrix" )
    labs <- as.character( labels[idx] )
    if( length( labs ) != length( pm.gid ) ) stop( "The number of elements in labels does not correspond to the number of genomes in the pan-matrix" )
  }
  if( show.labels ){
    text( Z[,x], Z[,y], labs, col=cols, ... )
  } else {
    points( Z[,x], Z[,y], pch=pch, col=cols, ... )
  }
}

plotLoadings <- function( pan.pca, x=1, y=2, show.labels=TRUE, col="black", pch=16, ... ){
  L <- pan.pca$Loadings
  xr <- range( L[,x] )
  yr <- range( L[,y] )
  args <- list(...)
  ii <- match( c("xlab","ylab"), names( args ) )
  if( is.na(ii[1]) ){
    xlab=paste( "PC", x, " (", round( 100*pan.pca$Evar[x] ), "%)", sep="" )
  } else {
    xlab <- args[[ii[1]]]
  }
  if( is.na(ii[2]) ){
    ylab=paste( "PC", y, "(", round( 100*pan.pca$Evar[y] ), "%)", sep="" )
  } else {
    ylab <- args[[ii[2]]]
  }
  
  plot( xr, c(0,0), type="l", col="gray", xlim=xr, ylim=yr, xlab=xlab, ylab=ylab )
  points( c(0,0), yr, type="l", col="gray" )
  
  if( show.labels ){
    text( L[,x], L[,y], rownames( L ), col=col, ... )
  } else {
    points( L[,x], L[,y], pch=pch, col=col, ... )
  }
}
