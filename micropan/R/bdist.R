#$Id: bdist.R 136 2013-07-20 08:33:59Z larssn $


bDist <- function( blast.files, e.value=1, verbose=TRUE ){

  if( verbose ) cat( "bDist:\n" )
  gids <- gregexpr( "GID[0-9]+", blast.files, extract=T )
  q.gid <- sapply( gids, function(x){ return( x[1] )} )
  h.gid <- sapply( gids, function(x){ return( x[2] )} )
  self.idx <- which( q.gid == h.gid )
  if( verbose ) cat( "...reading", length( self.idx ), "self alignments...\n" )
  q.max <- character( 0 )
  s.max <- numeric( 0 )
  qq <- character( 0 )
  hh <- character( 0 )
  dd <- numeric( 0 )
  for( i in 1:length( self.idx ) ){
    cat( "   ...reading file", blast.files[self.idx[i]], "\n" )
    tab <- readBlastTable( blast.files[self.idx[i]] )
    tab <- tab[which( tab$E.value <= e.value ),]
    tab <- tab[order( tab$Bit.score, decreasing=T ),]
    qh <- paste( tab$Query, tab$Hit )
    idx <- which( !duplicated( qh ) )
    q0 <- tab$Query[idx]
    h0 <- tab$Hit[idx]
    s0 <- tab$Bit.score[idx]
    
    q.max <- c( q.max, tapply( q0, q0, function(x){x[1]} ) ) # returns in sorted order
    s.max <- c( s.max, tapply( s0, q0, max ) ) # returns in sorted order

    d0 <- 1 - s0 / s.max[as.numeric( factor( q0, levels=q.max ) )]
 
    i1 <- which( q0 == h0 )
    i2 <- which( q0 != h0 )
    q0.s <- q0[i1]
    h0.s <- h0[i1]
    d0.s <- d0[i1]
    q0 <- q0[i2]
    h0 <- h0[i2]
    d0 <- d0[i2]
    
    qh <- paste( q0, h0 )
    hq <- paste( h0, q0 )
    nreci.idx <- which( !(qh %in% hq) ) #lacking reciprocal hits
    if( length( nreci.idx ) > 0 ){
      q0 <- c( q0, h0[nreci.idx] )
      h0 <- c( h0, q0[nreci.idx] )
      d0 <- c( d0, rep( 1, length( nreci.idx ) ) ) #distance 1 corresponds to no hit
      qh <- paste( q0, h0 )
    }
    ixx <- order( qh )
    q0 <- q0[ixx]
    h0 <- h0[ixx]
    d0 <- d0[ixx]  #must be sorted by qh here
    
    hq <- paste( h0, q0 )
    ixx.hq <- order( hq )
    d0 <- ( d0 + d0[ixx.hq] ) / 2
    idx.keep <- which( (1:length(d0)) < ixx.hq ) 
    qq <- c( qq, q0.s, q0[idx.keep] )
    hh <- c( hh, h0.s, h0[idx.keep] )
    dd <- c( dd, d0.s, d0[idx.keep] )
  }
  if( verbose ) cat( "...found BLAST results for", length( q.max ), "unique sequences...\n" )

  if( verbose ) cat( "...reading remaining alignments...\n" )
  u.gid <- unique( q.gid )
  cc <- length( self.idx )
  ctot <- length( blast.files )
  for( i in 1:(length( u.gid )-1) ){
    for( j in (i+1):length( u.gid ) ){
      idx <- which( (q.gid == u.gid[i]) & (h.gid == u.gid[j]) )
      if( length( idx ) == 1 ){
        tab <- readBlastTable( blast.files[idx] )
        tab <- tab[which( tab$E.value <= e.value ),]
        tab <- tab[order( tab$Bit.score, decreasing=T ),]
        qh <- paste( tab$Query, tab$Hit )
        idx <- which( !duplicated( qh ) )
        q0 <- tab$Query[idx]
        h0 <- tab$Hit[idx]
        s0 <- tab$Bit.score[idx]
      } else {
        stop( "Cannot find unique result file with genome", u.gid[i], "vs.", u.gid[j], ". Please supply a complete set of BLAST result files." )
      }
      
      idx <- which( (q.gid == u.gid[j]) & (h.gid == u.gid[i]) )
      if( length( idx ) == 1 ){
        tab <- readBlastTable( blast.files[idx] )
        tab <- tab[which( tab$E.value <= e.value ),]
        tab <- tab[order( tab$Bit.score, decreasing=T ),]
        qh <- paste( tab$Query, tab$Hit )
        idx <- which( !duplicated( qh ) )
        q0 <- c( q0, tab$Query[idx] )
        h0 <- c( h0, tab$Hit[idx] )
        s0 <- c( s0, tab$Bit.score[idx] )
      } else {
        stop( "Cannot find unique result file with genome", u.gid[j], "vs.", u.gid[i], ". Please supply a complete set of BLAST result files." )
      }
      
      d0 <- 1 - s0 / s.max[as.numeric( factor( q0, levels=q.max ) )]
  
      qh <- paste( q0, h0 )
      hq <- paste( h0, q0 )
      nreci.idx <- which( !(qh %in% hq) ) #lacking reciprocal hits
      if( length( nreci.idx ) > 0 ){
        q0 <- c( q0, h0[nreci.idx] )
        h0 <- c( h0, q0[nreci.idx] )
        d0 <- c( d0, rep( 1, length( nreci.idx ) ) ) #distance 1 corresponds to no hit
        qh <- paste( q0, h0 )
      }
      ixx <- order( qh )
      q0 <- q0[ixx]
      h0 <- h0[ixx]
      d0 <- d0[ixx]  #must be sorted by qh here
      
      hq <- paste( h0, q0 )
      ixx.hq <- order( hq )
      d0 <- ( d0 + d0[ixx.hq] ) / 2
      idx.keep <- which( (1:length(d0)) < ixx.hq )
      qq <- c( qq, q0[idx.keep] )
      hh <- c( hh, h0[idx.keep] )
      dd <- c( dd, d0[idx.keep] )
      
      cc <- cc + 2
      if( verbose ) cat( "  ...done with", cc, "out of", ctot, "files, have computed", length( dd ), "distances...\n" )
    }
  }
  ids <- which( qq == hh )
  utag <- unique( c( qq, hh ) )
  idn <- which( !(utag %in% qq[ids]) )
  if( length( idn ) > 0 ){
    if( verbose ) cat( "...adding some lacking self-alignments...\n" )
    qq <- c( qq, utag[idn] )
    hh <- c( hh, utag[idn] )
    dd <- c( dd, rep( 0, length( idn ) ) )
  } 
  dtab <- data.frame( Sequence.A=qq, Sequence.B=hh, Distance=dd, stringsAsFactors=F )
  return( dtab )
}
