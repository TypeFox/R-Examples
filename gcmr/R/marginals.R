## Marginals

## Marginal must have class marginal.gcmr and the following elements:
## - start(y, x, z, offset): compute initial estimates ignoring correlations
##      optional attributes lower and upper can be used
##      to specify box constrained parameters
## - dp(y, x, z, offset, lambda): evaluate [d,p]
## - q(p, x, z, offset, lambda): evaluate quantiles
## - npar(x, z): number of parameters
## - type: response type (integer or numeric)

# Gaussian
gaussian.marg <- function(link = "identity" ) {
    fm <- gaussian( substitute( link ) )
    ans <- list()
    ans$start <- function(y, x, z, offset) {
      if( !is.null(z) )
        offset <- list( as.vector(offset$mean), as.vector(offset$precision) )  
      eps <- sqrt(.Machine$double.eps)
      m <- glm.fit( x , y, offset=offset$mean, family=fm )
      sigma <- max( 10*eps, sd( residuals(m) ) )
      lambda <- c( coef(m), rep.int( 0, NCOL(z) ) )
      lambda[ NCOL(x)+1 ] <- ifelse( is.null(z), sigma, log(sigma) )

      if( is.null(z) ){
        names( lambda ) <- c( dimnames( as.matrix(x) )[[ 2L ]], "sigma" )
        attr( lambda, "lower" ) <- c( rep( -Inf, NCOL(x) ), eps )
      }
      else
        names( lambda ) <- c( paste("mean", dimnames( as.matrix(x) )[[ 2L ]], sep="."),
                             paste("dispersion", dimnames( as.matrix(z) )[[ 2L ]], sep=".") )
      lambda
    }
    ans$npar <- function(x, z) ifelse( !is.null(z), NCOL(x)+NCOL(z), NCOL(x)+1 )
    ans$dp <- function(y, x, z, offset, lambda) {
      nb <- length(lambda)
      mu <- fm$linkinv( x %*% lambda[ 1:NCOL(x) ] + offset$mean )
      if( is.null(z) )
        sd <- lambda[ nb ]
      else
        sd <- exp( z %*% lambda[ ( NCOL(x)+1 ):nb ] + offset$precision )
      cbind( dnorm( y , mu , sd ) , pnorm( y , mu , sd ) )
    }
    ans$q <- function(p, x, z, offset, lambda) {
      nb <- length(lambda)
      mu <- fm$linkinv( x %*% lambda[ 1:NCOL(x) ] + offset$mean )
      if( is.null(z) )
        sd <- lambda[ nb ]
      else
        sd <- exp( z %*% lambda[ ( NCOL(x)+1 ):nb ] + offset$precision )
      qnorm( p , mu , sd )
    }
    ans$fitted.val <- function(x, z, offset, lambda){
        fm$linkinv( x %*% lambda[ 1:NCOL(x) ] + offset$mean )
    }
    ans$type <- "numeric"
    class(ans) <- c( "marginal.gcmr")
    ans
}

# Binomial
# y is [#success,#failure]
binomial.marg <- function(link = "logit") {
    fm <- binomial( substitute( link ) )
    ans <- list()
    sizes  <- 1
    ans$start <- function(y, x, z, offset) {
        if(NCOL(y)==1) {
            y <- as.factor(y)
            y <- y!=levels(y)[1L]
            y <- cbind(y, 1-y)      
        } 
        sizes <<- y[,1]+y[,2]
        lambda <- coef( glm.fit( x, y, offset=offset$mean, family=fm ) )
        names(lambda) <- dimnames( as.matrix(x) )[[2L]]
        lambda
    }
    ans$npar <- function(x, z) NCOL(x)
    ans$dp <- function(y, x, z, offset, lambda) {
        if(NCOL(y)==1){
            y <- as.factor(y)
            y <- y!=levels(y)[1L]
            y <- cbind(y, 1-y)
        }
        mu <- fm$linkinv( x %*% lambda + offset$mean )
        cbind(dbinom( y[,1], sizes, mu ) ,
              pbinom( y[,1], sizes, mu ) )
    }
    ans$q <- function(p, x, z, offset, lambda) {
        mu <- fm$linkinv( x %*% lambda + offset$mean )
        q <- qbinom( p, sizes, mu )
        cbind( q, sizes-q )
    }
    ans$fitted.val <- function(x, z, offset, lambda){
        fm$linkinv( x %*% lambda + offset$mean )
    }
    ans$type <- "integer"
    class(ans) <- c( "marginal.gcmr")
    ans
}

# Poisson
poisson.marg <- function(link = "log") {
    fm <- poisson( substitute( link ) )
    ans <- list()
    ans$start <- function(y, x, z, offset) {
        lambda <- coef( glm.fit( x , y, offset=offset$mean, family=fm ) )
        names(lambda) <- dimnames( as.matrix(x) )[[ 2L ]]
        lambda
    }
    ans$npar <- function(x, z) NCOL(x)
    ans$dp <- function(y, x, z, offset, lambda) {
        mu <- fm$linkinv( x %*% lambda + offset$mean )
        cbind( dpois( y , mu ) , ppois( y , mu ) )
    }
    ans$q <- function(p, x, z, offset, lambda) {
        mu <- fm$linkinv( x %*% lambda + offset$mean )
        qpois( p , mu )
    }
    ans$fitted.val <- function(x, z, offset, lambda){
        fm$linkinv( x %*% lambda + offset$mean )
    }
    ans$type <- "integer"
    class(ans) <- c( "marginal.gcmr")
    ans
}


# Negative binomial
# var(y) = E(y) + k*E(y)^2 (k>0)
negbin.marg <- function(link = "log" ) {
  fm <- poisson( substitute( link ) )
  ans <- list()
  ans$start <- function(y, x, z, offset) {
      if( !is.null(z) )
        offset <- list( as.vector(offset$mean), as.vector(offset$precision) ) 
      eps <- sqrt(.Machine$double.eps)
      m <- glm.fit( x , y, offset=offset$mean, family=fm )
      mu <- fitted(m)
      kappa <- max( 10*eps , mean( ( (y-mu)^2-mu )/mu^2 ) )
      lambda <- c( coef(m), rep.int( 0, NCOL(z) ) )
      lambda[ NCOL(x)+1 ] <- ifelse( is.null(z), kappa, log(kappa) )
      if( is.null(z) ){
          names( lambda ) <- c( dimnames( as.matrix(x) )[[ 2L ]], "dispersion" )
          attr( lambda, "lower" ) <- c( rep(-Inf, NCOL(x) ), eps )
      }
      else
          names( lambda ) <- c( paste("mean", dimnames( as.matrix(x) )[[ 2L ]], sep="."),
                               paste("dispersion", dimnames( as.matrix(z) )[[ 2L ]], sep=".") )
      lambda
  }
  ans$npar <- function(x, z) ifelse( !is.null(z), NCOL(x)+NCOL(z), NCOL(x)+1 )
  ans$dp <- function(y, x, z, offset, lambda) {
      nb <- length(lambda)
      mu <- fm$linkinv( x %*% lambda[ 1:NCOL(x) ] + offset$mean )
      if( is.null(z) )
          size <- 1 / lambda[ nb ]
      else
          size <- 1 / exp( z %*% lambda[ ( NCOL(x)+1 ):nb ] + offset$precision )
      cbind( dnbinom( y, mu=mu, size=size) , pnbinom( y, mu=mu, size=size) )
  }
  ans$q <- function(p, x, z, offset, lambda) {
      nb <- length(lambda)
      mu <- fm$linkinv( x %*% lambda[ 1:NCOL(x) ] + offset$mean )
      if( is.null(z) )
          size <- 1 / lambda[ nb ]
      else
          size <- 1 / exp( z %*% lambda[ ( NCOL(x)+1 ):nb ] + offset$precision )
      qnbinom( p, mu=mu, size=size)
  }
  ans$fitted.val <- function(x, z, offset, lambda){
      fm$linkinv( x %*% lambda[ 1:NCOL(x) ] + offset$mean )
  }
  ans$type <- "integer"
  class(ans) <- c( "marginal.gcmr")
  ans
}

# next lines neeeded for back-compatibility with gcmr version 0.3
gs.marg <- gaussian.marg
bn.marg <- binomial.marg
ps.marg <- poisson.marg
nb.marg <- negbin.marg

# Weibull
weibull.marg <- function(link = "log"){
    fm <- Gamma( substitute( link ) ) # ;-)
    ans <- list()
    ans$start <- function(y, x, z, offset) {
      if( !is.null(z) )
        offset <- list( as.vector(offset$mean), as.vector(offset$precision) )  
      eps <- sqrt(.Machine$double.eps)
      m <- glm.fit(x , y, offset=offset$mean, family=fm)
      shape <- max( 10*eps, 1.2/sqrt( mean( log( y/fitted(m) )^2) ) )
      lambda <- c( coef(m), rep.int( 0, NCOL(z) ) )
      lambda[ NCOL(x)+1 ] <- ifelse( is.null(z), shape, log(shape) )

      if( is.null(z) ){
        names( lambda ) <- c( dimnames( as.matrix(x) )[[ 2L ]], "shape" )
        attr( lambda, "lower" ) <- c( rep( -Inf, NCOL(x) ), eps )
      }
      else
        names( lambda ) <- c( paste("scale", dimnames( as.matrix(x) )[[ 2L ]], sep="."),
                             paste("shape", dimnames( as.matrix(z) )[[ 2L ]], sep=".") )
      lambda
    }
    ans$npar <- function(x, z) ifelse( !is.null(z), NCOL(x)+NCOL(z), NCOL(x)+1 )
    ans$dp <- function(y, x, z, offset, lambda){
      nb <- length(lambda)
      scale <- fm$linkinv( x %*% lambda[ 1:NCOL(x) ] + offset$mean )
      if( is.null(z) )
        shape <- lambda[ nb ]
      else
        shape <- exp( z %*% lambda[ ( NCOL(x)+1 ):nb ] + offset$precision )
      cbind(dweibull(y, shape=shape, scale=scale) ,
            pweibull(y, shape=shape, scale=scale) )
    }
    ans$q <- function(p, x, z, offset, lambda){
      nb <- length(lambda)
      scale <- fm$linkinv( x %*% lambda[ 1:NCOL(x) ] + offset$mean )
      if( is.null(z) )
        shape <- lambda[ nb ]
      else
        shape <- exp( z %*% lambda[ ( NCOL(x)+1 ):nb ] + offset$precision )
      qweibull(p, shape=shape, scale=scale)
    }
    ans$fitted.val <- function(x, z, offset, lambda){
      fm$linkinv( x %*% lambda[ 1:NCOL(x) ] + offset$mean )
  }
    ans$type <- "numeric"
    class(ans) <- c( "marginal.gcmr")
    ans
} 
 
# Gamma
Gamma.marg <- function(link = "inverse"){
  fm <- Gamma( substitute( link ) )
  ans <- list()
  ans$start <- function(y, x, z, offset) {
    if( !is.null(z) )
      offset <- list( as.vector(offset$mean), as.vector(offset$precision) ) 
    eps <- sqrt(.Machine$double.eps)
    m <- glm.fit(x , y, offset=offset$mean, family=fm)
    disp <- sum( residuals(m, "pearson")^2 )/( NROW(y)-NCOL(x) )
    shape <- max( 10*eps, 1/disp )
    lambda <- c( coef(m), rep.int( 0, NCOL(z) ) )
    lambda[ NCOL(x)+1 ] <- ifelse( is.null(z), shape, log(shape) )

    if( is.null(z) ){
        names( lambda ) <- c( dimnames( as.matrix(x) )[[ 2L ]], "shape" )
        attr( lambda, "lower" ) <- c( rep( -Inf, NCOL(x) ), eps )
      }
      else
        names( lambda ) <- c( paste("mean", dimnames( as.matrix(x) )[[ 2L ]], sep="."),
                             paste("shape", dimnames( as.matrix(z) )[[ 2L ]], sep=".") )
    lambda
  }
  ans$npar <- function(x, z) ifelse( !is.null(z), NCOL(x)+NCOL(z), NCOL(x)+1 )
  ans$dp <- function(y, x, z, offset, lambda){
    nb <- length(lambda)
    mu <- fm$linkinv( x %*% lambda[ 1:NCOL(x) ] + offset$mean )
    if( is.null(z) )
      shape <- lambda[ nb ]
    else
      shape <- exp( z %*% lambda[ ( NCOL(x)+1 ):nb ] + offset$precision )
    cbind(dgamma(y, shape=shape, rate=shape/mu) ,
          pgamma(y, shape=shape, rate=shape/mu) )
  }
  ans$q <- function(p, x, z, offset, lambda){
    nb <- length(lambda)
    mu <- fm$linkinv( x %*% lambda[ 1:NCOL(x) ] + offset$mean )
    if( is.null(z) )
      shape <- lambda[ nb ]
    else
      shape <- exp( z %*% lambda[ ( NCOL(x)+1 ):nb ] + offset$precision )
    qgamma(p, shape=shape, rate=shape/mu)
  }
  ans$fitted.val <- function(x, z, offset, lambda){
      fm$linkinv( x %*% lambda[ 1:NCOL(x) ] + offset$mean )
  }
  ans$type <- "numeric"
  class(ans) <- c( "marginal.gcmr")
  ans
}

# beta
beta.marg <- function(link = "logit"){
    fm <- binomial( substitute( link ) ) # ;-)
    ans <- list()
    ans$start <- function(y, x, z, offset) {
        if( !is.null(z) )
            offset <- list( as.vector(offset$mean), as.vector(offset$precision) )   
        m <- betareg.fit(x=x, y=as.vector(y), z=z, offset=offset, link=link ) 
        lambda <- unlist( coef(m) )
        if( is.null(z) ){
            pos <- NCOL(x)+1
            lambda[pos] <- exp( lambda[pos] )
            names(lambda)[pos] <- "dispersion"
            attr(lambda, "lower") <- c( rep( -Inf, NCOL(x) ), sqrt(.Machine$double.eps) )
        }
    lambda
  }
  ans$npar <- function(x, z) ifelse(!is.null(z), NCOL(x)+NCOL(z), NCOL(x)+1)
  ans$dp <- function(y, x, z, offset, lambda) {
    nb <- length(lambda)
    mu <- fm$linkinv( x %*% lambda[ 1:NCOL(x) ] + offset$mean )
    if( is.null(z) )
      phi <- lambda[ nb ]
    else
      phi <- exp( z %*% lambda[ ( NCOL(x)+1 ):nb ] + offset$precision )
    shape1 <- phi*mu
    shape2 <- phi*(1-mu)
    cbind( dbeta(as.vector(y), shape1, shape2), pbeta(as.vector(y), shape1, shape2) )
  }
  ans$q <- function(p, x, z, offset, lambda) {
    nb <- length(lambda)
    mu <- fm$linkinv( x %*% lambda[1:NCOL(x)] + offset$mean )
    if( is.null(z) )
      phi <- lambda[ nb ]
    else
      phi <- exp( z %*% lambda[ ( NCOL(x)+1 ):nb ] + offset$precision )
    shape1 <- phi*mu
    shape2 <- phi*(1-mu)
    qbeta(p, shape1, shape2)
  }
  ans$fitted.val <- function(x, z, offset, lambda){
      fm$linkinv( x %*% lambda[ 1:NCOL(x) ] + offset$mean )
  }
  ans$type <- "numeric"
  class(ans) <- c("marginal.gcmr")
  ans
}



 
