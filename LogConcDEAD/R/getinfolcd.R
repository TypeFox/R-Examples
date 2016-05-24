### If we just have the values of the log-concave maximum likelihood estimator,
### this function gives more information

'getinfolcd' <- function(
                         x, ## data
                         y, ## logMLE
                         w=rep( 1/length(y), length(y) ),
                         chtol=10^-6,
                         MinSigma=NA,
                         NumberOfEvaluations=NA ){
  
  n <- nrow( x )
  d <- ncol( x )
### First find the points on the outside
  if (ncol(x) == 1) {
    chull=c(which.min(x), which.max(x) )}
  else {
    chull <- convhullnew(x)
  }
  
  ## points on the outside
  inhull <- c( unique( chull ) )
  
  xy <- cbind( x, y )
  xymin <- cbind( x[  inhull, ],  min( y ) - 1 )
  allpoints <- rbind( xy, xymin )

  ## First of all, find the relevant bunique and betaunique
  chopt <- paste( "C", format( chtol, scientific=FALSE ), sep="" )
  smoothch <- convhullnew( allpoints, options=chopt )
  smoothch <- smoothch[ apply( smoothch <= n, 1, sum ) == d + 1 , , drop=FALSE ]
  nbbetaunique <- nrow( smoothch )
  bunique <- array( dim=c( nbbetaunique, d ) )
  betaunique <- rep( 0, nbbetaunique )
  dropme <- rep( FALSE, nbbetaunique )
  for( j in 1:nbbetaunique ) {
    A <-  t( x[ smoothch[ j, -1],] ) - x[ smoothch[j,1],]
    detA <- abs(det( A ))
    if( detA > 10^-8 ) {
      alpha <- x[ smoothch[ j, 1 ], ]
      z <- y[ smoothch[ j, -1 ] ] - y[ smoothch[ j, 1 ] ]
      bunique[ j, ] <- solve( t( A ), z )
      betaunique[ j ] <-  alpha %*% bunique[ j, ] - y[ smoothch[ j, 1 ] ]
    } else {
      dropme[ j ] <- TRUE
    }
  }
  if( sum( dropme ) )  {
    betaunique <- betaunique[!dropme]
    bunique <- bunique[!dropme,]
  }
  nfree <- length( unique( c( smoothch[!dropme,] ) ) )
  
  ## Secondly, do it for the triangulation
  chopt <- paste( "Qt C", format( chtol, scientific=FALSE ), sep="" )
  triang <- convhullnew( allpoints, options=chopt )
  triang<- triang[ apply( triang <= n, 1, sum ) == d + 1 , , drop=FALSE ]
  nrows <- nrow( triang )
  verts <- array(dim=c(nrows,d,d))
  vertsoffset <- array(dim=c(nrows,d))
  A <- array( dim=c( nrows, d, d ) )
  alpha <- array( dim=c( nrows, d ) )
  detA <- rep( 0, nrows )
  dropme <- rep( FALSE, nrows )
  b <- array( dim=c( nrows, d ) )
  beta <- rep( 0, nrows )
  for (j in 1:nrows) {
    A[j, ,] <- t( x[ triang[ j, -1],] ) - x[ triang[j,1],] 
    detA[ j ] <- abs( det( as.matrix( A[j, ,] ) ) )
    if( detA[ j ] > 10^-8 ) {
      alpha[j, ] <- x[ triang[ j, 1 ], ]
      z <- y[ triang[j,-1]] - y[triang[j,1]]
      b[j,] <- solve( t( A[j,,] ), z )
      beta[ j ] <- alpha[j,] %*% b[j,] - y[triang[j,1]]
      verts[j,,] <- solve( A[j,,] )
      vertsoffset[j,] <- verts[j,,]%*% alpha[j,]
    } else {
      dropme[ j ] <- TRUE
    }
  }
  
  rownames(b) <- NULL
  ## deal with the flat ones
  if( sum( dropme ) )  {
    triang <- triang[!dropme,]
    beta <- beta[!dropme]
    b <- b[!dropme,]
    A <- A[!dropme,,]
    alpha <- alpha[!dropme,]
    verts <- verts[!dropme,,]
    vertsoffset <- vertsoffset[!dropme,]
    detA <- detA[!dropme]
  }

  ## Get the outer stuff
  if(ncol(x)==1) {
    midpoint = mean(x)
    outnorm <- matrix(c(-1,1),ncol=1)
    outdist <- abs( range( x ) - midpoint )
  }
  else {
    outnorm <- matrix( NA, nrow=nrow(chull), ncol=d )
    outdist <- rep( NA, nrow( chull ) )
    midpoint <- apply( x, 2, mean) ##find a midpoint
    for (i in 1:nrow(chull)) {
      tmp <- t( apply( x[ chull[i,] ,] , 1, "-", midpoint ) )
      tmp <- solve( tmp, rep( 1, d ) )
      norm <- sqrt( sum( tmp^2 ) )
      outnorm[ i, ] <- tmp/norm
      outdist[ i ] <- 1/norm
    } 
  }

  r <- list(
            x = x,
            w = w,
            logMLE = y,
            NumberOfEvaluations = NumberOfEvaluations,
            MinSigma = MinSigma,
            b = b,
            beta = beta,
            triang = triang,
            verts = verts,
            vertsoffset = vertsoffset,
            chull = chull,
            outnorm = outnorm,
            outdist = outdist,
            midpoint = midpoint,
            A = A,
            alpha =  alpha,
            detA = detA,
            bunique = bunique,
            betaunique = betaunique,
            nfree = nfree )
  class(r) <- "LogConcDEAD"
  return(r)
}
