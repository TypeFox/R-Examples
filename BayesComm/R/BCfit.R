BCfit <-
function(y, X, covlist, R, z, mu, updateR, iters, thin = 1, burn = 0, priW = c(nrow(z) + 2 * ncol(z), 2 * ncol(z)), verbose = 0) {
  
  spnames <- colnames(y)
  y <- t(y)         # NOTE: y is transposed relative to the input!
  nsp <- dim(y)[1]  # number of species
  n <- dim(y)[2]    # number of sites
  iR <- solve(R)    # Inverse correlation matrix
  e <- z - mu       # component of z that can't be explained by mu
  
  nsamp <- (iters - burn) %/% thin
  
  # Dave asks: Should mu be in the output as well? It wasn't before, but that 
  #    may have been an oversight.
  # Nick: It's not needed since it's just X %*% B, we could
  #   store it for convenience but I doubt it's worth the memory used
  output <- list(
    R = array(NA, dim = c(nsamp, ((nsp * nsp - nsp) / 2))),
    B = NULL, 
    z = array(NA, dim = c(nsamp, n, nsp)), 
    burn = burn, 
    thin = thin
  )
  
  for (i in 1:nsp) {
    temp <- matrix(NA, nsamp, length(covlist[[i]]))
    colnames(temp) <- colnames(X)[covlist[[i]]]
    output$B[[spnames[i]]] <- temp
  }
  rm(temp)
  
  nam <- rep(NA, n * n)
  for (i in 1:nsp) {
    for (j in 1:nsp) {
      nam[(i - 1) * nsp + j] <- paste(spnames[i], "_", spnames[j], sep="")
    }
  }
  
  colnames(output$R) <- nam[which(upper.tri(diag(nsp)))]
  dimnames(output$z)[[3]] <- spnames
  
  # start sampler
  rec <- 0 # counter for record number after burn-in and thinning
  start <- Sys.time()
  for (iter in 1:iters) {
    
    # get truncation values and sample z
    trunc <- find_trunc(mu, y)
    e <- sample_e(e, trunc, iR)
    z <- mu + e
    
    # sample mu and calculate e
    mulis<- sample_mu(z, X, covlist)
    mu <- mulis[[1]]
    e <- z - mu
    
    # sample R
    if (updateR) {
      R <- sample_R(z - mu, priW)
      iR <- chol2inv(chol(R))
    }
    
    if(verbose == 2){
      message(iter)
    }
    if(verbose > 0 & iter == burn){
      message("burn-in complete")
    }
    
    # record parameters
    if (iter %% thin == 0 & iter > burn) {
      if(verbose == 1){
        message(iter)
      }
      rec <- rec + 1
      output$R[rec, ] <- R[upper.tri(R)]
      output$z[rec, , ] <- z
      for (i in 1:nsp) {
        output$B[[i]][rec, ] <- mulis[[2]][[i]] 
      }
    }
  }  # sampler
  
  output
}
