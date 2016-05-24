
## =============================================================================
## jacobian.full  : generates a full jacobian matrix by numerical differencing
## =============================================================================

jacobian.full<- function(y, func, dy=NULL, time=0,
     parms=NULL, pert=1e-8, ...) {

  if (!is.numeric(y))
    stop("y-values should be numeric")

  N    <- length(y)
  refy <- y
  ifelse (is.null(dy), refdy <- unlist( func(time,y,parms,...))[1:N], refdy <- dy)

  if (! is.numeric(refdy))
    stop("dy-values should either be NULL or numeric")
  if (length(refdy) != N)
    stop("function should return at least one value for each y")

  ynames  <- attr(y,"names")
 
# Perturb the state variables one by one
  delt   <- perturb(y,pert)

  ny <-length(y)
  jacob  <- matrix(nrow=ny,ncol=ny,data=0)
  for (j in 1:ny) {
    y[j] <- y[j]+delt[j]

    # recalculate model rate of change
    dy  <-  unlist( func(time,y,parms,...))[1:N]

    # impact of the current state var on rate of change of all state vars
    jacob [,j] <- (dy-refdy)/delt[j]

    y[j] <- refy[j]   # restore
  }
  colnames (jacob) <- ynames
           
  return(jacob)

}
