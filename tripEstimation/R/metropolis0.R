"metropolis0" <-
function(model,
                       iters=1000,thin=10,
                       start.x=NULL,start.z=NULL) {
  
  ## Initialize x,z
  x0 <- start.x
  z0 <- start.z
  if(is.null(x0)) x0 <- model$start.x
  if(is.null(z0)) z0 <- model$start.z

  if (is.null(x0)) stop("no starting points for X")
  if (is.null(z0)) stop("no starting points for Z")
  
  ## Drop dimnames for speed
  dimnames(x0) <- NULL
  dimnames(z0) <- NULL
  
  ## Extract model components
  proposal.x <- model$proposal.x
  proposal.z <- model$proposal.z
  mask.x <- model$mask.x
  mask.z <- model$mask.z
  logp.position <- model$logp.position
  logp.behavioural <- model$logp.behavioural
  fixed.x <- model$fixed.x
  
  ## Number of locations
  n <- nrow(x0)
  
  ## Allocate storage for the chain
  ch.x <- array(0,c(dim(x0),iters))
  ch.z <- array(0,c(n-1,2,iters))
  
  ## Contribution to logp from the initial x.
  logp.posn0 <- logp.position(x0)
  mask.x0 <- mask.x(x0)
  mask.z0 <- mask.z(z0)
  for(k1 in 1:iters) {
    for(k2 in 1:thin) {

      ## Propose all x at once, and calculate mask and contribution to
      ## the log posterior
      x1 <- proposal.x(x0)
      x1[fixed.x,] <- x0[fixed.x,]
      mask.x1 <- mask.x(x1)
      logp.posn1 <- logp.position(x1)


      ## Update x
      ## In each case we compute full contribution (positional +
      ## behavourial) to the log posterior for current and proposed
      ## points, and apply the MH rule. If the proposal is accepted,
      ## we update both x and the cached positional contribution to
      ## the log posterior.
      
      ## Accept/reject first x
      if(mask.x1[1] || !mask.x0[1]) {
        logp0 <- logp.posn0[1]+logp.behavioural(1,x0[1,1:2,drop=FALSE],z0[1,,drop=FALSE],x0[2,1:2,drop=FALSE])
        logp1 <- logp.posn1[1]+logp.behavioural(1,x1[1,1:2,drop=FALSE],z0[1,,drop=FALSE],x1[2,1:2,drop=FALSE])
        if(logp1-logp0 > log(runif(1))) {
          x0[1,] <- x1[1,]
          logp.posn0[1] <- logp.posn1[1]
          mask.x0[1] <- mask.x1[1]
        }
      }
      
      ## Accept/reject last x
      if(mask.x1[n] || !mask.x0[n]) {
        logp0 <- logp.posn0[n]+logp.behavioural(n-1,x0[n-1,1:2,drop=FALSE],z0[n-1,,drop=FALSE],x0[n,1:2,drop=FALSE])
        logp1 <- logp.posn1[n]+logp.behavioural(n-1,x1[n-1,1:2,drop=FALSE],z0[n-1,,drop=FALSE],x1[n,1:2,drop=FALSE])
        if(logp1-logp0 > log(runif(1))) {
          x0[n,] <- x1[n,]
          logp.posn0[n] <- logp.posn1[n]
          mask.x0[n] <- mask.x1[n]
        }
      }

      ## Red/Black update for interior x
      for(rb in 2:3) {
        is <- seq(rb,n-1,by=2)
        is <- is[mask.x1[is] | !mask.x0[is]]
        logp0 <- (logp.posn0[is]+
                  logp.behavioural(is-1,x0[is-1,1:2,drop=FALSE],z0[is-1,,drop=FALSE],x0[is,1:2,drop=FALSE])+
                  logp.behavioural(is,x0[is,1:2,drop=FALSE],z0[is,,drop=FALSE],x0[is+1,1:2,drop=FALSE]))
        logp1 <- (logp.posn1[is]+
                  logp.behavioural(is-1,x1[is-1,1:2,drop=FALSE],z0[is-1,,drop=FALSE],x1[is,1:2,drop=FALSE])+
                  logp.behavioural(is,x1[is,1:2,drop=FALSE],z0[is,,drop=FALSE],x1[is+1,1:2,drop=FALSE]))
        ## MH rule - compute indices of the accepted points.
        accept <- is[logp1-logp0 > log(runif(length(is)))]
        x0[accept,] <- x1[accept,]
        logp.posn0[accept] <- logp.posn1[accept]
        mask.x0[accept] <- mask.x1[accept]
      }
      
      ## Update z
      ## Here we need only consider the behavioural contributions to
      ## the log posterior (the position contributions are constant
      ## and would cancel), and so we can update all the z in parallel.
      z1 <- proposal.z(z0)
      mask.z1 <- mask.z(z1)
      is <- (1:(n-1))[mask.z1 | !mask.z0]
      logp0 <- logp.behavioural(is,x0[is,1:2,drop=FALSE],z0[is,],x0[is+1,1:2,drop=FALSE])
      logp1 <- logp.behavioural(is,x0[is,1:2,drop=FALSE],z1[is,],x0[is+1,1:2,drop=FALSE])
      ## MH rule - compute indices of the accepted points.
      accept <- is[(logp1-logp0 > log(runif(length(is))))]
      z0[accept,] <- z1[accept,]
      mask.z0[accept] <- mask.z1[accept]
    }
    ## Store the current state
    ch.x[,,k1] <- x0
    ch.z[,,k1] <- z0
  }
  
  ## ensure masks in result reflect fixed positions
  mask.x0[fixed.x] <- TRUE
  
  if (all(mask.x0)) print("All X's pass the mask") else {print("X's do not pass mask:");print( which(!mask.x0))}
  if (all(mask.z0)) print("All Z's pass the mask") else {print("Z's do not pass mask:"); print( which(!mask.z0))}
  list(model=model,
       x=ch.x,z=ch.z,
       mask.x=mask.x0,mask.z=mask.z0,
       last.x=x0,last.z=z0)
}

