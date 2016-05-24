dwpt.boot <- function(y, wf, J=log(length(y),2)-1, p=1e-04, frac=1) {

  N <- length(y)
  if(N/2^J != trunc(N/2^J)) 
    stop("Sample size is not divisible by 2^J")
  
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
  y.basis <- as.logical(ortho.basis(portmanteau.test(y.dwpt, p, type="other")))

  ## Taken from my 2D bootstrapping methodology
  resample.dwpt <- y.dwpt
  for(i in 1:length(y.basis)) {
    m <- length(y.dwpt[[i]])
    if(y.basis[i])
      resample.dwpt[[i]] <- sample(y.dwpt[[i]], replace=TRUE)
    else
      resample.dwpt[[i]] <- rep(NA, m)
  }
  idwpt(resample.dwpt, y.basis)
}

