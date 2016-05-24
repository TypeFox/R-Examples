# Author: D.R. Bickel
# Modications: P. Poncet
grenander <-
function(x,        # sample (the data)
         bw = NULL,# fraction of the observations to be considered
         k,        # parameter 'k'
         p,        # parameter 'p' (power);  if p==Inf, 'mlv.venter' is used
         ...)
{
############################
# Grenander's mode estimator
############################
  
  if (p == Inf) {
    cat("argument 'p' is infinite. Venter's mode estimator is used")
    return(venter(x = x, bw = bw, k = k, ...)) 
  }

  ny <- length(x)    

  if (missing(k) & !is.null(bw)) {
    if (bw <= 0 | bw > 1) stop("argument 'bw' must belong to (0, 1]")
    k <- ceiling(bw*ny) - 1
  } else if (missing(k) & is.null(bw)) {
    k <- ceiling(ny/2) - 1
  }
  
  if (k < 0 | k >= ny) stop("argument 'k' must belong to [0, length('x'))") 

  y <- sort(x)

  inf <- y[1:(ny-k)]
  sup <- y[(k+1):ny]
  diff <- sup - inf
  tot <- inf + sup
  if(any(diff==0)) {
    warning("limiting value of Grenander mode used") #! ??
    M <- mean(ifelse(diff==0, tot, NA), na.rm = TRUE)/2
  } else {
    b <- sum(tot/diff^p)/2
    a <- sum(1/diff^p)
    if(is.finite(b/a)) {
      M <- b/a
    } else {
      stop("function 'grenander' failed. Argument 'p' may be too large")
    }
  }
  
  ## Output
  return(M)
}

#Grenander <- grenander
