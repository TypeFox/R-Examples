#$Id: rarefaction.R 166 2014-07-13 21:41:15Z khliland $


rarefaction <- function( pan.matrix, n.perm=1 ){
  pan.matrix[which( pan.matrix > 0, arr.ind=T )] <- 1
  ng <- dim( pan.matrix )[1]
  nmat <- matrix( 0, nrow=ng, ncol=n.perm )
  cm <- apply( pan.matrix, 2, cumsum )
  nmat[,1] <- rowSums( cm>0 )
  if( n.perm > 1 ){
    cat( "permuting:\n" )
    for( i in 2:n.perm ){
      cm <- apply( pan.matrix[sample( ng ),], 2, cumsum )
      nmat[,i] <- rowSums( cm>0 )
      cat( "." )
      if( (i/100)==round(i/100) ) cat( "\n" )
    }
    cat( "\n" )
  }
  rownames( nmat ) <- paste( 1:ng, "genomes" )
  colnames( nmat ) <- paste( "Permutation", 1:n.perm )
  class( nmat ) <- c( "Rarefac", "matrix" )
  return( nmat )
}



plot.Rarefac <- function( x, type="b", pch=16, xlab="Genomes", ylab="Number of unique gene clusters", ... ){
  Rarefac <- x
  plot( 1:dim( Rarefac )[1], rowMeans( Rarefac ), type=type, pch=pch, xlab=xlab, ylab=ylab, ... )
}

summary.Rarefac <- function( object, ... ){
  Rarefac <- object
  cat( "For", 1, "genome we observe on average", round( mean( Rarefac[1,] ) ), "unique gene clusters\n" )
  for( i in 2:dim( Rarefac )[1] ){
    cat( "For", i, "genomes we observe on average", round( mean( Rarefac[i,] ) ), "unique gene clusters\n" )
  }
}

str.Rarefac <- function( object, ... ){
  Rarefac <- object
  cat( "Rarefaction over", dim( Rarefac )[1], "genomes\n" )
}