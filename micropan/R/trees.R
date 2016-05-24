#$Id: trees.R 170 2014-07-23 12:57:06Z larssn $


panTree <- function( pan.matrix, dist.FUN=distManhattan, nboot=0, linkage="average", ... ){
  P <- ncol( pan.matrix )
  distFUN <- match.fun( dist.FUN )
  htree <- hclust( distFUN( pan.matrix, ... ), method=linkage )
  
  nbranch=NULL
  if( nboot > 0 ){
    cat( "bootstrapping" )
    signatur <- clusterSignature( htree$merge )
    nbranch <- rep( 0, length( signatur ) )
    names( nbranch ) <- signatur
    for( i in 1:nboot ){
      cat( "." )
      ht <- hclust( distFUN( pan.matrix[,sample( (1:P), P, replace=T )], ... ), method=linkage )
      nbranch <- nbranch + as.numeric( signatur %in% clusterSignature( ht$merge ) )
    }
    cat( "\n" )
  }
  pantree <- list( Htree=htree, Nboot=nboot, Nbranch=nbranch, Dist.FUN=as.character( substitute( dist.FUN ) ) )
  class( pantree ) <- c( "Pantree", "list" )
  return( pantree )
}


plot.Pantree <- function( x, leaf.lab=NULL, col="black", xlab="", main="", cex=1, show.boot=TRUE, ... ){
  # x is a Pantree
  # labels...
  if( is.null( leaf.lab ) ){
    labs <- x$Htree$labels
  } else {
    if( is.null( names( leaf.lab ) ) ) stop( "Each element in leaf.lab must be named by its GID.tag" )
    pm.gid <- x$Htree$labels
    gid <- names( leaf.lab )
    idx <- match( pm.gid, gid )
    if( sum( is.na( idx ) ) > 0 ) stop( "GID.tag names does not match GID.tags in the tree" )
    labs <- as.character( leaf.lab[idx] )
    if( length( labs ) != length( pm.gid ) ) stop( "The number of elements in leaf.lab does not correspond to the number of genomes in the tree" )
  }
  # colors...
  if( length( col )==1 ){
    col.lab <- rep( col, length.out=length( labs ) )
  } else {
    if( is.null( names( col ) ) ) stop( "Each element in col must be named by its GID.tag" )
    pm.gid <- x$Htree$labels
    gid <- names( col )
    idx <- match( pm.gid, gid )
    if( sum( is.na( idx ) ) > 0 ) stop( "GID.tag names does not match GID.tags in the tree" )
    col.lab <- col[idx]
    if( length( col ) != length( pm.gid ) ) stop( "The number of elements in col does not correspond to the number of genomes in the tree" )
  }
  x$Htree$labels <- labs
  dendro <- dendrapply( as.dendrogram( x$Htree ), setLeafAttributes, labs, col.lab, cex )
  
  cpar <- par()$mar
  maxnc <- max( nchar( labs ) )
  mar4 <- 2 + cex*(maxnc*0.8*par("cin")[1]*2.54)
  par( mar=c(5,1,1,mar4) )
  plot( dendro, horiz=T, xlab=xlab, main=main, ... )
  if( (x$Nboot > 0) & show.boot ){
    lab.nod <- as.character( round( 100*x$Nbranch/x$Nboot )/100 )
    bp <- branchPos( x$Htree$merge, x$Htree$order )
    for( i in 1:length( x$Htree$height ) ){
      if( !(x$Htree$merge[i,1]<0 & x$Htree$merge[i,2]<0) ){
        text( x$Htree$height[i], bp[i], lab.nod[i], cex=0.75, col="red4", pos=4, offset=0.1 )
      }
    }
  }
  par( mar=cpar )  
}


summary.Pantree <- function( object, ... ){
  # object is a Pantree
  labs <- object$Htree$labels
  cat( "Pangenome tree for ", length( labs ), " genomes (", paste( labs[1:3], collapse="," ), "...)\n", sep="" )
  cat( "Distances computed by", object$Dist.FUN, "and containing", object$Nboot, "bootstrap samples\n")
}


str.Pantree <- function( object, ... ){
  # object is a Pantree
  labs <- object$Htree$labels
  cat( "Pangenome tree for", length( labs ), "genomes\n" )
}


setLeafAttributes <- function( node, lab.names, lab.col, lab.cex ){
  if( is.leaf( node ) ){
    attr( node, "nodePar" ) <- list( lab.cex=lab.cex, pch=NA, lab.col=lab.col[which( lab.names == attr( node, "label" ) )] )
  }
  return( node )
}


clusterSignature <- function( mergeMatrix ){
  N <- nrow( mergeMatrix )
  signature <- character( N )
  for( i in 1:N ){
    if( mergeMatrix[i,1] < 0 ){
      left <- as.character( -1*mergeMatrix[i,1] )
    } else {
      left <- gsub( ";", ",", signature[mergeMatrix[i,1]] )
    }
    if( mergeMatrix[i,2] < 0 ){
      right <- as.character( -1*mergeMatrix[i,2] )
    } else {
      right <- gsub( ";", ",", signature[mergeMatrix[i,2]] )
    }
    left <- paste( sort( unlist( strsplit( left, split="," ) ) ), collapse="," )
    right <- paste( sort( unlist( strsplit( right, split="," ) ) ), collapse="," )
    zig <- sort( c( left, right ) )
    signature[i] <- paste( zig[1], zig[2], sep=";" )
  }
  return( signature )
}


branchPos <- function( mergeMatrix, ordering ){
  N <- nrow( mergeMatrix )
  branchp <- numeric( N )
  for( i in 1:N ){
    if( mergeMatrix[i,1] < 0 ){
      left <- which( ordering == -1*mergeMatrix[i,1] )
    } else {
      left <- branchp[mergeMatrix[i,1]]
    }
    if( mergeMatrix[i,2] < 0 ){
      right <- which( ordering == -1*mergeMatrix[i,2] )
    } else {
      right <- branchp[mergeMatrix[i,2]]
    }
    branchp[i] <- mean( c( left, right ) )
  }
  return( branchp )
}
