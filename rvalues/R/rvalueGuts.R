
## Auxiliary functions
rvalueGuts <- function(dat, alpha.grid, V, vfun, hypers, smooth) {
  ########################################################
  # dat           nunits x 2 [specific to family]
  # alpha.grid    ngrid
  # V             nunits x ngrid
  # vfun          function tail prob
  # hypers        length 2; hyper-parameters [specific to familiy in vfun]
  ########################################################
  
  nunits <- nrow(dat)
  ngrid <- length(alpha.grid)
  cc <- numeric(ngrid)
  for( j in 1:ngrid )
  {
      cc[j] <- quantile(V[,j], prob= 1 - alpha.grid[j], names = FALSE, type = 1)
  }
  
  ## smooth and functionalize
  if(smooth=="none") {
      ccfun <- approxfun(alpha.grid, cc, yleft = 1, yright = 0)
  }
  else {
      cc2 <- supsmu(alpha.grid, cc, bass=smooth )
      ccfun <- approxfun(c(0,cc2$x,1), c(1,cc2$y,0))
  }
  ### Think of ccfun as the lambda_{\alpha} function
  
  #dfun <- function(alpha, unitdata, hypers)
  #{
 #   dd <- ccfun(alpha) - vfun(alpha, unitdata, hypers)
  #  dd
  #}
  ## march through units finding rvalue by uniroot
  ### What if there are multiple roots? (This may be possible if ccfun is not monotonic)
  #rvals <- numeric(nunits)
  #for( i in 1:nunits ) {
  #   rvals[i] <- uniroot(dfun,interval=c(0,1),unitdata=dat[i,],hypers=hypers, tol=1e-10)$root
  #}
  #rvals <- mroot(dfun, lower=rep(0,nunits), upper=rep(1,nunits), unitdata=dat, 
  #               hypers=hypers, tol=1e-6)$root
  rvals <- VVcut(V, cc, nunits, ngrid, alpha.grid)
  ## maybe go back to root finding?
  #rvals <- alpha.grid[tmp]
  
  ans <- list()
  ans$rvals <- rvals
  ans$lamfun <- cc
  ans$smoothlamfun <- ccfun
  return(ans)
}


vfun.pg <- function( alpha, unitdata, hypers ) {
  ## Poisson Gamma posterior upper tail probability
  aa <- hypers[1]
  bb <- hypers[2]  ## hyper params of Gamma prior
  #x <- unitdata[,1]
  #et <- unitdata[,2]	
  x <- unitdata[1]
  et <- unitdata[2]
  tAlpha <- qgamma(alpha, shape=aa, rate=bb, lower.tail=FALSE )
  p <- pgamma(tAlpha, shape=(aa+x), rate=(bb+et), lower.tail=FALSE )
  p
}

vfun.bb <- function( alpha, unitdata, hypers ) {
  aa <- hypers[1]
  bb <- hypers[2]
  if(is.matrix(unitdata)) {
     x <- unitdata[,1]
     n <- unitdata[,2]
  }
  else {
     x <- unitdata[1]
     n <- unitdata[2]
  }
  tAlpha <- qbeta(alpha, shape1=aa, shape2=bb, lower.tail=FALSE )
  p <- pbeta(tAlpha, shape1=(aa+x), shape2=(bb+n-x), lower.tail=FALSE )
  p
}

vfun.nn <- function( alpha, unitdata, hypers ) {
  ## ignors hypers...already converted to theta ~ Normal(0,1)
  tPM <- unitdata[,1]/(1+unitdata[,2])
  tSD <- sqrt( unitdata[,2]/(1+unitdata[,2]) )
  tAlpha <- qnorm(alpha, lower.tail=FALSE )
  z <- (tAlpha  - tPM)/tSD
  p <- pnorm(z, lower.tail=FALSE )
  p
}

vfun.gg <- function(alpha, unitdata, hypers) {
  aa <- hypers[1]
  bb <- hypers[2]  ## hyper params of Gamma prior
  x <- unitdata[,1]
  ss <- unitdata[,2]
	
  tAlpha.inv <- qgamma(alpha, shape=aa, scale=bb )
  p <- pgamma(tAlpha.inv, shape = aa + ss, rate = x + bb)
  p
}


