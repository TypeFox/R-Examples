is.dummy <- function( v )
{
  ( sum( (v!=0) & (v!=1) ) == 0 ) &
  ( sum( v!=1 ) != 0 )
}

get.scale <- function( X )
{
  s.scale <- apply( X, 2, sd )
  s.scale[1] <- 1.0
  for( i in 2:dim(X)[2] )
    if( is.dummy( X[,i] ) )
      s.scale[i] = 1.0	 
  return( s.scale )
}

zic <- function( formula, data, a0, b0, c0, d0, e0, f0, n.burnin, n.mcmc, n.thin, tune = 1.0, scale = TRUE ) 
{
  # unsorted data matrices 
  mdl <- model.frame( formula, data )
  y <- model.response( mdl )
  X <- model.matrix( formula, mdl )

  # sort matrices	 
  idx <- sort( y, index.return = TRUE )$ix
  y <- y[idx]
  X <- X[idx,]

  if( scale )
    {
      s.scale <- get.scale( X )
      X <- scale( X, center = FALSE, scale = s.scale )
    }

  # call C++
  output <- .Call( "zic_sample", 
                   y, X, 
	           a0, b0, -9, -9, -9, -9, -9, e0, f0, 
	           c0, d0, -9, -9, -9, -9, -9,
	           FALSE, n.burnin, n.mcmc, n.thin, tune, package = "zic" )
 
  output$alpha <- mcmc( output$alpha )
  output$beta <- mcmc( output$beta )	
  output$gamma <- mcmc( output$gamma )	
  output$delta <- mcmc( output$delta )	
  output$sigma2 <- mcmc( output$sigma2 )	

  varnames(output$alpha) <- list( "alpha" )
  varnames(output$beta) <- colnames(X)[2:dim(X)[2]]
  varnames(output$gamma) <- list( "gamma" )
  varnames(output$delta) <- colnames(X)[2:dim(X)[2]]
  varnames(output$sigma2) <- list( "sigma2" )

  if( scale )
    {
      output$s.scale <- s.scale[2:length(s.scale)]
    }
 
  return( output )
}

zic.svs <- function( formula, data, a0, g0.beta, h0.beta, nu0.beta, r0.beta, s0.beta, e0, f0, 
                                    c0, g0.delta, h0.delta, nu0.delta, r0.delta, s0.delta, 
                                    n.burnin, n.mcmc, n.thin, tune = 1.0, scale = TRUE )
{
  # unsorted data matrices 
  mdl <- model.frame( formula, data )
  y <- model.response( mdl )
  X <- model.matrix( formula, mdl )

  # sort matrices	 
  idx <- sort( y, index.return = TRUE )$ix
  y <- y[idx]
  X <- X[idx,]

  if( scale )
    {
      s.scale <- get.scale( X )
      X <- scale( X, center = FALSE, scale = s.scale )
    }
  
  # call C++	      
  output <- .Call( "zic_sample", 
                   y, X, 
		   a0, -9, g0.beta, h0.beta, nu0.beta, r0.beta, s0.beta, e0, f0, 
		   c0, -9, g0.delta, h0.delta, nu0.delta, r0.delta, s0.delta,
                   TRUE, n.burnin, n.mcmc, n.thin, tune, package = "zic" )

  output$alpha <- mcmc( output$alpha )	   
  output$beta <- mcmc( output$beta )	
  output$gamma <- mcmc( output$gamma )	
  output$delta <- mcmc( output$delta )	
  output$sigma2 <- mcmc( output$sigma2 )	
  output$I.beta <- mcmc( output$I.beta )	
  output$I.delta <- mcmc( output$I.delta )	
  output$omega.beta <- mcmc( output$omega.beta )
  output$omega.delta <- mcmc( output$omega.delta )

  varnames(output$alpha) <- list( "alpha" )
  varnames(output$beta) <- colnames(X)[2:dim(X)[2]]
  varnames(output$gamma) <- list( "gamma" )
  varnames(output$delta) <- colnames(X)[2:dim(X)[2]]
  varnames(output$sigma2) <- list( "sigma2" )
  varnames(output$I.beta) <- colnames(X)[2:dim(X)[2]]
  varnames(output$I.delta) <- colnames(X)[2:dim(X)[2]]
  varnames(output$omega.beta) <- list( "omega.beta" )
  varnames(output$omega.delta) <- list( "omega.delta" )

  if( scale )
    {
      output$s.scale <- s.scale[2:length(s.scale)]
    }
  
  return( output )
}


