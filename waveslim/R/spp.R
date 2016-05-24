spp.mle <- function(y, wf, J=log(length(y),2)-1, p=0.01, frac=1)
{
  sppLL <- function(x, y) {
    delta <- x[1]
    fG <- x[2]

    ## cat("Parameters are: d =", delta, ", and f =", fG, fill=TRUE)
    
    y.dwpt <- y[[1]]
    y.basis <- y[[2]]
    n <- y[[3]]
    J <- y[[4]]
    
    ## Establish the limits of integration for the band-pass variances
    a <- unlist(apply(matrix(2^(1:J)-1), 1, seq, from=0, by=1)) / 
      2^(rep(1:J, 2^(1:J))) / 2
    b <- unlist(apply(matrix(2^(1:J)), 1, seq, from=1, by=1)) / 
      2^(rep(1:J, 2^(1:J))) / 2
    
    ## Define some useful parameters for the wavelet packet tree
                                        # n <- length(y)
    length.jn <- n / rep(2^(1:J), 2^(1:J))
    scale.jn <- rep(2^(1:J+1), 2^(1:J))
    
    ## Initialize various parameters for the reduced LL
    Basis <- (1:length(y.basis))[y.basis]
    bp.var <- numeric(length(Basis))
    delta.n <- 100

    ## Compute the band-pass variances according to \delta and f_G
    omega.diag <- NULL
    for(i in 1:sum(y.basis)) {
      jn <- Basis[i]
      bp.var[i] <- bandpass.spp(a[jn], b[jn], delta, fG)
      omega.diag <- c(omega.diag,
                      scale.jn[jn] * rep(bp.var[i], length.jn[jn]))
    }
    
    ## Compute reduced log-likelihood 
    rLL <- n * log(1/n * sum(y.dwpt^2 / omega.diag, na.rm=TRUE)) +
      sum(length.jn[y.basis] * log(scale.jn[y.basis] * bp.var))
    rLL
  }
  
  n <- length(y)
  x0 <- numeric(2)

  ## Perform discrete wavelet packet transform (DWPT) on Y
  y.dwpt <- dwpt(y, wf, n.levels=J)
  n <- length(y)
  if(frac < 1) {
    for(i in 1:length(y.dwpt)) {
      vec <- y.dwpt[[i]]
      ni <- length(vec)
      j <- rep(1:J, 2^(1:J))[i]
      vec[trunc(frac * n/2^j):ni] <- NA
      y.dwpt[[i]] <- vec
    }
  }
  y.basis <- as.logical(ortho.basis(portmanteau.test(y.dwpt, p)))
  y.dwpt <- as.matrix(unlist(y.dwpt[y.basis]))

  ## Compute initial estimate of the Gegenbauer frequency
  y.per <- per(y - mean(y))
  x0[2] <- (0:(n/2)/n)[max(y.per) == y.per]

  ## Compute initial estimate of the fractional difference parameter
  muJ <- (unlist(apply(matrix(2^(1:J)-1), 1, seq, from=0, by=1)) / 
          2^(rep(1:J, 2^(1:J))) + 
          unlist(apply(matrix(2^(1:J)), 1, seq, from=1, by=1)) /  
          2^(rep(1:J, 2^(1:J)))) / 4
  y.modwpt <- modwpt(y, wf=wf, n.levels=J)
  y.varJ <- rep(2^(1:J), 2^(1:J)) *
    unlist(lapply(y.modwpt,
                  FUN=function(x)sum(x*x,na.rm=TRUE)/length(x[!is.na(x)])))
  x0[1] <- min(-0.5 * lsfit(log(abs(muJ[y.basis] - x0[2])), 
                            log(y.varJ[y.basis]))$coef[2], 0.49)

  cat(paste("Initial parameters are: delta =", round(x0[1],4), 
            "freqG =", round(x0[2],4), "\n"))
  result <- optim(par=x0, fn=sppLL, method="L-BFGS-B",
                  lower=c(0.001,0.001), upper=c(0.499,0.499),
                  control=list(trace=0, fnscale=2),
                  y=list(y.dwpt, y.basis, n, J))
  return(result)
}

spp2.mle <- function(y, wf, J=log(length(y),2)-1, p=0.01,
                     dyadic=TRUE, frac=1)
{

  spp2LL <- function(x, y) {
    d1 <- x[1]
    f1 <- x[2]
    d2 <- x[3]
    f2 <- x[4]
    ## cat("Parameters are: d1 =", round(d1,6), ", and f1 =", round(f1,6),
    ##     ", d2 =", round(d2,6), ", and f2 =", round(f2,6), fill=TRUE)

    y.dwpt <- y[[1]]
    y.basis <- y[[2]]
    n <- y[[3]]
    J <- y[[4]]

    ## Establish the limits of integration for the band-pass variances
    a <- unlist(apply(matrix(2^(1:J)-1), 1, seq, from=0, by=1)) / 
      2^(rep(1:J, 2^(1:J))) / 2
    b <- unlist(apply(matrix(2^(1:J)), 1, seq, from=1, by=1)) / 
      2^(rep(1:J, 2^(1:J))) / 2
    
    ## Define some useful parameters for the wavelet packet tree
    length.jn <- n / rep(2^(1:J), 2^(1:J))
    scale.jn <- rep(2^(1:J+1), 2^(1:J))
    
    ## Initialize various parameters for the reduced LL
    Basis <- (1:length(y.basis))[y.basis]
    bp.var <- numeric(length(Basis))
    delta.n <- 100

    ## Compute the band-pass variances according to \delta and f_G
    omega.diag <- NULL
    for(i in 1:sum(y.basis)) {
      jn <- Basis[i]
      bp.var[i] <- bandpass.spp2(a[jn], b[jn], d1, f1, d2, f2)
      omega.diag <- c(omega.diag,
                      scale.jn[jn] * rep(bp.var[i], length.jn[jn]))
    }
    
    ## Compute reduced log-likelihood 
    n * log(1/n * sum(y.dwpt^2 / omega.diag, na.rm=TRUE)) +
      sum(length.jn[y.basis] * log(scale.jn[y.basis] * bp.var), na.rm=TRUE)
  }

  n <- length(y)
  x0 <- numeric(4)

  ## Perform discrete wavelet packet transform (DWPT) on Y
  y.dwpt <- dwpt(y, wf, n.levels=J)
  if(!dyadic) {
    for(i in 1:length(y.dwpt)) {
      vec <- y.dwpt[[i]]
      ni <- length(vec)
      j <- rep(1:J, 2^(1:J))[i]
      vec[trunc(frac * n/2^j):ni] <- NA
      y.dwpt[[i]] <- vec
    }
  }
  y.basis <- as.logical(ortho.basis(portmanteau.test(y.dwpt, p, type="other")))
  y.dwpt <- as.vector(unlist(y.dwpt[y.basis]))

  ## Compute initial estimate of the Gegenbauer frequencies
  if(dyadic)
    y.per <- per(y - mean(y))
  else
    y.per <- per(y[1:(frac*n)] - mean(y[1:(frac*n)]))
  freq.y <- (0:(frac*n %/% 2))/(frac*n)
  x0[2] <- freq.y[max(y.per) == y.per]
  x0[4] <- freq.y[max(y.per[freq.y > x0[2] + freq.y[10] | 
    freq.y < x0[2] - freq.y[10]]) == y.per]
  if(x0[2] > x0[4]) {
    xx <- x0[2]
    x0[2] <- x0[4]
    x0[4] <- xx
    rm(xx)
  }

  ## Compute initial estimate of the fractional difference parameters
  muJ <- (unlist(apply(matrix(2^(1:J)-1), 1, seq, from=0, by=1)) / 
      2^(rep(1:J, 2^(1:J))) + 
      unlist(apply(matrix(2^(1:J)), 1, seq, from=1, by=1)) /  
      2^(rep(1:J, 2^(1:J)))) / 4
  y.modwpt <- modwpt(y, wf=wf, n.levels=J)
  y.varJ <- rep(2^(1:J), 2^(1:J)) *
    unlist(lapply(y.modwpt,
                  FUN = function(x) sum(x*x,na.rm=TRUE)/length(x[!is.na(x)])))
  x0.mid <- (x0[2] + x0[4]) / 2
  muJ <- muJ[y.basis]
  y.varJ <- y.varJ[y.basis]
  x0[1] <- min(-0.5 * lsfit(log(abs(muJ[muJ < x0.mid] - x0[2])), 
                            log(y.varJ[muJ < x0.mid]))$coef[2], 0.49)
  x0[3] <- min(-0.5 * lsfit(log(abs(muJ[muJ > x0.mid] - x0[4])), 
                            log(y.varJ[muJ > x0.mid]))$coef[2], 0.49)

  cat(paste("Initial parameters: d1 = ", round(x0[1],4), 
            ", f1 = ", round(x0[2],4), ", d2 = ", round(x0[3],4), 
            ", f2 = ", round(x0[4],4), sep=""), fill=TRUE)
  result <- optim(par=x0, fn=spp2LL, method="L-BFGS-B",
                  lower=rep(0.001,4), upper=rep(0.499,4),
                  control=list(trace=1, fnscale=2),
                  y=list(y.dwpt, y.basis, n, J))
  return(result)
}
