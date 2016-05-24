#$Id: eclust.R 179 2014-08-23 20:52:52Z khliland $


bClust <- function( dist.table, linkage="single", threshold=1.0 ){
  cat( "bClust:\n" )
  linknum <- grep( linkage, c( "single", "average", "complete" ) )
  
  dt <- dist.table[which( dist.table$Distance < threshold ),]
  utag <- sort( unique( c( dt$Sequence.A, dt$Sequence.B ) ) ) # Important to sort here!
    
  cat( "...constructing graph with", length( utag ), "sequences (nodes) and", dim( dt )[1], "distances (edges)\n" )
  M <- matrix( as.numeric( factor( c( dt$Sequence.A, dt$Sequence.B ), levels=utag ) ), ncol=2, byrow=F )
  g <- graph.edgelist( M, directed=F )
  cls <- clusters( g )
  cat( "...found", cls$no, "single linkage clusters\n" )
  clustering <- cls$membership + 1
  names( clustering ) <- utag
  
  if( linknum > 1 ){
    ucls <- unique( cls$membership )
    incomplete <- sapply( 1:cls$no, function( j ){
      v <- which( cls$membership == ucls[j] )
      degg <- degree( g, v )
      return( min( degg ) < (length( degg ) + 1) )
    })
    cat( "...found", sum( incomplete ), "incomplete clusters, splitting:\n")
    clustering <- clustering * 1000
    inc <- which( incomplete )     #the incomplete clusters
    if( length( inc ) > 0 ){
      idx.inc <- which( cls$membership %in% inc )
      memnum <- cls$membership[idx.inc]
      memtag <- utag[idx.inc] #is also sorted since idx.inc and utag are sorted according to utag
      idi <- which( (dt$Sequence.A %in% utag[idx.inc]) & (dt$Sequence.B %in% utag[idx.inc]) )
      aa <- dt$Sequence.A[idi]
      bb <- dt$Sequence.B[idi]
      dd <- dt$Distance[idi]
      for( i in 1:length( inc ) ){
        idx <- which( memnum == inc[i] )
        dmat <- matrix( 1, nrow=length( idx ), ncol=length( idx ) )
        rownames( dmat ) <- memtag[idx]
        colnames( dmat ) <- memtag[idx]
        idd <- which( (aa %in% memtag[idx]) | (bb %in% memtag[idx]) )
        a <- as.numeric( factor( aa[idd], levels=memtag[idx] ) )
        b <- as.numeric( factor( bb[idd], levels=memtag[idx] ) )
        dmat[matrix( c(a,b), ncol=2, byrow=F )] <- dd[idd]
        dmat[matrix( c(b,a), ncol=2, byrow=F )] <- dd[idd]
        if( linknum == 2 ){
          clst <- hclust( as.dist( dmat ), method="average" )
        } else {
          clst <- hclust( as.dist( dmat ), method="complete" )
        }
        clustering[idx.inc[idx]] <- clustering[idx.inc[idx]] + cutree( clst, h=threshold )
        cat( "." )
        if( (i/100)==round(i/100) )cat( "\n" )
        aa <- aa[-idd]
        bb <- bb[-idd]
        dd <- dd[-idd]
      }
    }
    cat( "\n" )
  }
  cat( "...ended with", length( unique( clustering ) ), "clusters, largest cluster has", max( table( clustering ) ), "members\n" )
  clustering <- sort( clustering )
  return( clustering )
}



isOrtholog <- function( clustering, dist.table ){
  aa <- dist.table$Sequence.A
  bb <- dist.table$Sequence.B
  dd <- dist.table$Distance
  uclst <- unique( clustering )
  tags <- names( clustering )
  is.ortholog <- rep( F, length( clustering ) )
  names( is.ortholog ) <- tags
  for( i in 1:length( uclst ) ){
    idx <- clustering == uclst[i]
    seqz <- tags[idx]
    ns <- length( seqz )
    idd <- (aa %in% seqz) & (bb %in% seqz) 
    gidz <- sapply( gregexpr( "GID[0-9]+", seqz, extract=T ), function(x){return(x[1])} )
    if( max( table( gidz ) ) > 1 ){
      dmat <- matrix( 1, nrow=ns, ncol=ns )
      a <- as.numeric( factor( aa[idd], levels=seqz ) )
      b <- as.numeric( factor( bb[idd], levels=seqz ) )
      dmat[matrix( c(a,b), ncol=2, byrow=F )] <- dd[idd]
      dmat[matrix( c(b,a), ncol=2, byrow=F )] <- dd[idd]
      ixx <- order( rowSums( dmat ) )
      ixd <-  !duplicated( gidz[ixx] ) 
      is.ortholog[idx][ixx[ixd]] <- T
    } else {
      is.ortholog[idx] <- T
    }
#    aa <- aa[-idd]
#    bb <- bb[-idd]
#    dd <- dd[-idd]
    cat( "." )
    if( (i/100)==round(i/100) )cat( "." )
  }
  cat( "\n" )
  return( is.ortholog )
}
   
    
    
