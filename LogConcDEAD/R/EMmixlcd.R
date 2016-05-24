## Main function to compute EM mixture

'EMmixlcd' <- function( x, k = 2, y, props, epsratio=10^-6, max.iter=50,
epstheta=10^-8, verbose=-1 ) 
{
  x <- as.matrix( x )
  n <- nrow(x)
  d <- ncol(x)

  if( d==1 ) {
    stop("1D data, Please see package logcondens")
  }

  if( missing( y ) || missing( props ) ) {
    ## By default, a Gaussian hierarchical clustering mixture is used
    ## as a starting point - allow variances to vary
    require( mclust )
    highclust <- hc( modelName="VVV", x )
    class <- c( hclass( highclust, k ) )
    props <- rep( 0, k )
    y <- matrix( 0, nrow=n, ncol=k )
    for( i in 1:k ) {
      props[ i ] <- sum( class==i ) / n
      ss <- x[ class==i, ]
      y[ , i ] <- dmvnorm( x, mean=apply( ss,2, mean), sigma=var( ss ), log=TRUE )
    }
  }

  oldloglik <- logliklcdmix( y, props )
  hcy <- y

  ## Get the weights and remove duplicates
  xdup <- duplicated( x )
  tmp <- getweights( x )
  x <- tmp$x
  w <- tmp$w
  y <- matrix( 0, nrow=nrow( x ), ncol = k )
  for( i in 1:k ) {
    y[ ,i ] <- hcy[ !xdup, i]
  }

  ## update likelihood
  niter <- 1
  while( niter==1 || (abs(( newloglik - oldloglik ) / oldloglik ) > epsratio  &&
        niter <= max.iter )) {
    pif <- t( t( exp( y ) ) * props )
    theta <- pif / apply( pif, 1, sum )
    props <- apply( theta, 2, sum ) / sum( theta )

    for ( i in 1:k ) {
      whichx <- ( theta[ , i ]/sum( theta[ , i ] ) ) > epstheta / n
      tmp <- mlelcd( x[ whichx, ], w = w[ whichx ] * theta[ whichx, i ] / sum( w[ whichx ] * theta[ whichx, i ] ))
      y[ whichx, i ] <- tmp$logMLE
      y[ !whichx, i ] <- -Inf
    }

    props <- apply( theta, 2, sum ) / sum( theta )

    ## update likelihood
    if (niter > 1) oldloglik <- newloglik
    newloglik <- logliklcdmix( y, props )
    if(verbose > 0) cat( paste( "Number of iterations:", niter," Oldloglik:", oldloglik, "Newloglik",
                         newloglik, "\n" ) )
    niter <- niter + 1
  }

  if(niter > max.iter && verbose >= 0) warning( "maximum number of iterations reached" )
  r <- list( x = x,
            logf = y,
            props = props,
            niter = niter,
            lcdloglik = newloglik
            )
  class( r ) <- "lcdmix"
  return( r )
}
