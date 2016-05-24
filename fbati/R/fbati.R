## Get the actual coding of the markers
xcode <- function( m0, m1, model ) {
  if( model==ADDITIVE )
    return( as.integer(m0==2) + as.integer(m1==2) )
  if( model==DOMINANT )
    return( as.integer(m0==2 | m1==2) )
  if( model==RECESSIVE )
    return( as.integer(m0==2 & m1==2) )
  stop( paste( "xcode: model (", model, ") not understood.", sep="" ) )
}

## Get the 3 codings of the markers
xcodes <- function(model) {
  Xc <- c();
  if( model==ADDITIVE ) {
    Xc <- c(0,1,2);
  }else if( model==DOMINANT ){
    Xc <- c(0,1,1);
  }else if( model==RECESSIVE ){
    Xc <- c(0,0,1);
  }else{
    stop( paste("xcodes: model (",model,") is not understood.", sep="") );
  }

  return( Xc );
}

#############################
## calculating group means ##
#############################
groupMean <- function( group, x ) {
  ugroup <- unique(group);
  xbar <- rep( 0, length(group) );
  for( u in ugroup ){
    wh <- group==u;
    umean <- mean( x[wh] );
    xbar[wh] <- umean;
  }

  return(xbar);
}

## Only uses the first affected (traces to dataComputeGroupG C function),
##  so this is true for the LL method, DR could be different
fbati.calc <- function( data, ## what was passed to datamatrix
                        m0pos=7, m1pos=m0pos+1,
                        groupsG,  ## part of the result of datamatrix
                        affectedIndex, ## the other part of the result
                        envCol=m1pos+1,
                        model, iter=1000,
                        debug=FALSE ) {
  ## actually - we really first need to get rid of the data that isn't needed,
  ##  _before_ the groups are computed! (unaffected are also marked by group...)
  ## keep (1) children that (2) are affected

  ## The fix for when we have missing markers in the affected child...
  #keep <- groupsG$groups!=0 & data$AffectionStatus==2
  keep <- groupsG$groups!=0 & data$AffectionStatus==2 & data[,m0pos]!=0 & data[,m1pos]!=0
  keep <- intersect( which(keep), affectedIndex ) ## Modification so only uses first affected
  data <- data[keep,]
  groupsG <- groupsG[keep,]

  if( debug ) {
    cat( "*** data[keep,] ***\n" )
    print( data )
    cat( "*** groupsG ***\n" )
    print( groupsG )
  }

  ## get X(g), then the means
  x <- xcode( data[,m0pos], data[,m1pos], model )
  xbar <- groupMean( groupsG$groups, x )
  zbar <- groupMean( groupsG$groups, data[,envCol] )

  xmxbar <- x-xbar
  zmzbar <- data[,envCol]-zbar

  return( fbati2( xmxbar, zmzbar, groupsG$groups, iter, debug ) )
}

fbati2 <- function( xmxbar, zmzbar, group, iter=1000, debug=FALSE ){
  if( length(xmxbar)==1 ) return(1)

  o <- order( group );  ## NEEDS TO BE SORTED FIRST

  if( debug ) {
    print( data.frame( xmxbar=xmxbar[o], zmzbar=zmzbar[o], group=group[o] ) )
  }

  #print( data.frame( xmxbar=xmxbar, zmzbar=zmzbar, group=group ) )

  ## then call the C function for speed.
  pvalue <- as.double(0.0);
  res = .C( "fbati",
      pvalue,

      as.integer(length(group)),
      as.double(xmxbar[o]),
      as.double(zmzbar[o]),
      as.integer(group[o]),

      as.integer(iter),

      DUP=TRUE); #DUP=FALSE );
  pvalue = res[[1]] ## 03/21/2014 <-- Make sure that this works!

  #print( pvalue )

  numInf <- sum( xmxbar!=0 & zmzbar!=0 )
  strataSum <- strataSum( xmxbar, zmzbar, group )

  #return( pvalue );
  return( list( pvalue=pvalue, numInf=numInf, strataSum=strataSum ) )
}

strataSum <- function( xmxbar, zmzbar, group ) {
  o <- order( group )
  xmxbar <- xmxbar[o]
  zmzbar <- zmzbar[o]
  group <- group[o]

  ssum <- NULL
  for( g in unique(group) ) {
    wh <- group==g
    if( any(xmxbar[wh]!=0 & zmzbar[wh]!=0) ) {
      numInf=sum(xmxbar[wh]!=0 & zmzbar[wh]!=0)
      names(numInf) <- g; #paste( "G", g, sep="" )

      if( is.null(ssum) ) {
        ssum <- numInf
      }else{
        ssum <- c( ssum, numInf )
      }
    }
  }

  res <- t(data.frame(ssum))
  row.names(res) <- NULL
  return( res )
}

strataNameDehash <- function( number ) {
  str <- as.character("                                                  ")  ## 50 spaces...
  res <- .C( "pG_group_dehash", as.integer(number), str )
  ##print( res )
  return( res[[2]] )
}
