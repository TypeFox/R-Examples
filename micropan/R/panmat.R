
#$Id: panmat.R 167 2014-07-15 12:33:45Z larssn $

panMatrix <- function( clustering ){
  gids <- sapply( gregexpr( "GID[0-9]+", names( clustering ), extract=T ), function(x){x[1]} )
  ugids <- sort( unique( gids ) )
  ngids <- length( ugids )
  uclst <- sort( unique( clustering ) )
  nclst <- length( uclst )
  pan.matrix <- matrix( 0, nrow=ngids, ncol=nclst )
  rownames( pan.matrix ) <- ugids
  colnames( pan.matrix ) <- paste( "Cluster", uclst, sep="_" )

  for( i in 1:ngids ){
    idx <- which( gids == ugids[i] )
    clst <- clustering[idx]
    tab <- table( clst )
    idd <- as.numeric( names( tab ) )
    ixx <- which( uclst %in% idd )
    pan.matrix[i,ixx] <- tab
  }
  attr( pan.matrix, "clustering" ) <- clustering
  class( pan.matrix ) <- c( "Panmat", "matrix" )
  return( pan.matrix )
}


plot.Panmat <- function( x, col="black", xlab="Number of genomes", ylab="Number of clusters", ... ){
  # x is a Panmat
  x[which( x > 0, arr.ind=T )] <- 1
  levs <- 1:dim( x )[1]
  y <- table( factor( colSums( x ), levels=levs ) )
  barplot( y, col=col, border=col, names.arg=levs, xlab=xlab, ylab=ylab, ... )
}

summary.Panmat <- function( object, ... ){
  # object is a Panmat
  object[which( object > 0, arr.ind=T )] <- 1
  levs <- 1:dim( object )[1]
  y <- table( factor( colSums( object ), levels=levs ) )
  for( i in 1:length( levs ) ){
    cat( y[i], "clusters found in", levs[i], "genomes\n" )
  }
}

str.Panmat <- function( object, ... ){
  # object is a Panmat
  dd <- dim( object )
  cat( "Pan-matrix with", dd[2], "gene families over", dd[1], "genomes\n" )
}
