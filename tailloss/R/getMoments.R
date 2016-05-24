getMoments <- function (x, theta = 0, cap = Inf, maxmom = 10) {

  stopifnot(x >= 0,
            is.scalar(theta), theta >= 0,
            is.scalar(cap), cap > 0,
            is.scalar(maxmom), maxmom >= 1)
  maxmom <- floor(maxmom) # integer moments only

  if (theta == 0) {
  	
  	 if (!is.infinite(cap))
  	 	x <- pmin(cap, x)

    robj <- outer(log(x), 1:maxmom)
    
    robj <- exp(robj)

  } else { # theta > 0

    alpha <- 1 / theta^2
    beta <- alpha / x

    if (is.infinite(cap)) {

      robj <- outer(beta, 1:maxmom, function(b, k) {
        lgamma(alpha + k) - (k * log(b) + lgamma(alpha))
      })
      
        robj <- exp(robj)

    } else { # cap < Inf
    
      u <- cap
      p <- pgamma(u, shape = alpha, rate = beta)
      
      robj <- outer(beta, 1:maxmom, function(b, k) {
        igamma(alpha + k, b * u, log = TRUE) - (k * log(b) + lgamma(alpha))
      })
      
      robj <- exp(robj) + exp(outer(log(1 - p), log(u) * (1:maxmom), '+'))
    
    }
    
  }
  

  robj
}