##' initialize data structure for bootstrapping
##' 
##' this will create an array for the climate data and corresponding
##' matrix for the tree-ring data, holding the randomly resampled data
##' for bootstrapping
##' @param u climate data (matrix with parameters in columns and years in rows)
##' @param g vector with tree-ring data
##' @param n times for resampling
##' @param boot one of c("stationary", "std", "exact", "none")
##' @importFrom np b.star
##' @importFrom boot tsboot
##' @return a list
##' @keywords internal
init_boot_data <- function(u, g, n, boot) {
  m <- length(g)
  k <- dim(u)[2]
  
  if (boot == "none") {
    out_u <- u
    out_g <- g
  }
  
  if (boot == "stationary") {
    if (length(g) < 18) {
      stop("For stationary bootstrap, win_size has to be at least 18. Consider adapting win_size or changing the bootstrap procedure.")
    }
    ## following Politis and Romano 1994, and especially the automatic
    ## calculation of the optimal block length
    out_u <- array(dim = c(m, k, n))
    out_g <- matrix(nrow = m, ncol = n)
    ind <- 1:m

    b_star <- round(b.star(g)[1])
    if (b_star == 0) {
      b_star <- 1
      warning("Optimal block length was estimated to be smaller than 1 and set to 1 by treeclim. Usage of standard, non-stationary bootstrap is recommended in that case.")
    } 
    inx <- tsboot(tseries = ind, statistic = function(x) x,
                R = n, sim = "geom", l = b_star, orig.t = TRUE)$t

    for (i in 1:n) {
      out_u[, , i] <- u[inx[i, ], ]
      out_g[, i] <- g[inx[i, ]]
    }
  }
  
  if (boot %in% c("std", "dendroclim")) {
    out_u <- array(dim = c(m, k, n))
    out_g <- matrix(nrow = m, ncol = n)
    for (i in 1:n) {
      .sample <- sample(1:m, m, replace = TRUE)
      out_u[, , i] <- u[.sample, ]
      out_g[, i] <- g[.sample]
    }
  }
  
  if (boot == "exact") {
    ## gaussian simulation of tree-ring data with circulant embedding; this code 
    ## is adapted from Dave Meko's seascorr MATLAB source
    
    ## climate data remains unchanged
    out_u <- u
    
    ## subtract series mean
    trm <- g - mean(g)
    
    ## store original variance and sd
    varo <- var(g)
    xsd <- sd(g)
    
    ## taper with cosine bell
    trt <- spec.taper(trm, 0.05)
    
    ## rescale so that mean is exactly zero
    trt <- trt - mean(trt)
    
    ## find power of 2 greater than double length of data
    ll <- length(trt)
    pow2 <- 2^c(1:12)
    pow <- pow2[which(pow2 > 2 * ll)[1]]
    
    ## pad data with 0 to length of next power of 2
    padlen <- pow - ll
    trp <- c(trt, rep(0, padlen))
    
    ## scale variance of trp to original variance of u
    sdx <- sd(g)
    sdp <- sd(trp)
    meanp <- mean(trp)
    trtemp <- trp - meanp
    trtemp <- trtemp * (sdx / sdp)
    trp <- trtemp + meanp
    
    ## compute discrete fourier transform on tapered and padded series
    z <- fft(trp)
    
    ## compute periodogram; this differs from MATLAB version according
    ## to a scaling factor of ~ 2
    Pyy  <- Re(z * Conj(z))/padlen
    
    ## sample Gaussian noise and compute mu for 1000 simulation runs
    M <- length(trp)/2
    Z <- matrix(rnorm(500 * length(trp) * 2), ncol = 500)
    k <- t(1:length(trp))
    i1 <- 2 * k - 1
    i2 <- 2 * k
    
    term1 <- matrix(complex(pow, Z[i1,], Z[i2,]), ncol = 500)
    term2 <- sqrt(Pyy/length(trp))
    term2 <- matrix(rep(term2, 500), ncol = 500)
    mu <- term1 * term2
    
    V <- fft(mu)
    
    Vr <- Re(V)
    Vi <- Im(V)
    D1 <- Vr[1:M,]
    D2 <- Vi[1:M,]
    D <- cbind(D1, D2)
    D <- D[1:ll, 1:1000]
    
    ## scale to original mean and variance
    Dmean <- matrix(rep(colMeans(D), each = ll), ncol = 1000)
    Dsd <- matrix(rep(apply(D, 2, sd), each = ll), ncol = 1000)
    Dz <- (D - Dmean) / Dsd
    out_g <- Dz * matrix(xsd, nrow = ll, ncol = 1000) + mean(g)
    
  }
  
  list(climate = out_u,
       chrono = out_g)
}
