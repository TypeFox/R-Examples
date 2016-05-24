#$Id: domainclust.R 101 2013-05-24 16:43:19Z larssn $

dClust <- function( hmmer.table ){
  cat( "dClust:\n" )
  cat( "...hmmer.table contains", length( unique( hmmer.table$Query ) ), "proteins...\n" )
  cat( "...with hits against", length( unique( hmmer.table$Hit ) ), "HMMs...\n" )
  hmmer.table <- hmmer.table[order( hmmer.table$Start ),]
  dseq <- sort( unlist( tapply( hmmer.table$Hit, hmmer.table$Query, function( x ){ paste( x, collapse="," ) } ) ) )
  seq <- names( dseq )
  dsc <- as.numeric( factor( as.vector( dseq ) ) )
  names( dsc ) <- seq
  cat( "...ended with", length( unique( dsc ) ), "clusters, largest cluster has", max( table( dsc ) ), "members\n" )
  attr( dsc, "cluster.info" ) <- unique( as.vector( dseq ) )
  return( dsc )
}


hmmerCleanOverlap <- function( hmmer.table ){
  qt <- table( hmmer.table$Query )
  cat( "There are", length( qt ), "proteins in this hmmer.table...\n" )
  hmmer.table <- hmmer.table[order( hmmer.table$Start ),]
  idx <- which( qt > 1 )
  multi <- names( qt )[idx]
  cat( "There are", length( multi ), "proteins with multiple hits, resolving overlaps:\n")
  keep <- rep( T, dim( hmmer.table )[1] )
  for( i in 1:length( multi ) ){
    idx <- which( hmmer.table$Query == multi[i] )
    keep[idx] <- nonoverlap( hmmer.table[idx,] )
    if( (i/100) == round(i/100) ) cat( "." )
  }
  cat( "\n" )
  return( hmmer.table[keep,] )
}


nonoverlap <- function( hmmer.table ){
  nh <- dim( hmmer.table )[1]
  keep <- rep( T, nh )
  ht <- hmmer.table[keep,]
  dif <- ht$Start[2:nh] - ht$Stop[1:(nh-1)]
  ido <- which( dif <= 0 )
  while( (nh > 1) & (length( ido ) > 0) ){
    idx <- unique( c( ido, ido+1 ) )
    idd <- which( ht$Evalue[idx] == max( ht$Evalue[idx] ) )
    map <- which( keep )
    keep[map[idx[idd[1]]]] <- F
    ht <- hmmer.table[keep,]
    nh <- dim( ht )[1]
    if( nh > 1 ){
      dif <- ht$Start[2:nh] - ht$Stop[1:(nh-1)]
      ido <- which( dif <= 0 )
    }
  }
  return( keep )
}












