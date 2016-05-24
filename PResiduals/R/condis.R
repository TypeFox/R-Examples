condis2 <- function(A, B) {
    ## The input is the two numeric vectors to be compared

    n <- length(A)

    ## Create a matrix of data counts
    P <- matrix(0L, n, n)

    ## Create list of indexes for instertions
    P[cbind(rank(A),rank(B))] <- 1L

    nrow <- dim(P)[1]
    ncol <- dim(P)[2]
    if(nrow<2 | ncol<2) return(0)
  
    N <- sum(P)  # sample size

    ## P is now a probability distribution.
    P <- P/N
  
    ## P0 is the expected distribution under independence.
    prow <- P %*% rep(1,ncol)
    pcol <- rep(1, nrow) %*% P
    P0 <- prow %*% pcol
  
    Pbar <- P0bar <- matrix(,nrow,ncol)
  
    ## For each position (i,j), calculate the probabilities for the 4 sections:
    ## aa1[j] = probability for the section (<i, <j)
    ## aa2[j] = probability for the section (<i, >j)
    ## bb1[j] = probability for the section (>i, <j)
    ## bb2[j] = probability for the section (>i, >j)
    aa <- aa0 <- rep(0,ncol)
    bb <- bb0 <- pcol
    for(i in seq_len(nrow)) {
        if(i > 1) {
            aa0 <- aa0 + P0[i-1,]
            aa <- aa + P[i-1,]
        }
        
        bb0 <- bb0 - P0[i,]
        aa01 <- c(0, cumsum(aa0[-ncol]))
        aa02 <- sum(aa0) - cumsum(aa0)
        bb01 <- c(0, cumsum(bb0[-ncol]))
        bb02 <- sum(bb0) - cumsum(bb0)
        P0bar[i,] <- aa01 + bb02 - aa02 - bb01

        bb <- bb - P[i,]
        aa1 <- c(0, cumsum(aa[-ncol]))
        aa2 <- sum(aa) - cumsum(aa)
        bb1 <- c(0, cumsum(bb[-ncol]))
        bb2 <- sum(bb) - cumsum(bb)
        Pbar[i,] <- aa1 + bb2 - aa2 - bb1
    }
  
    u <- Pbar %*% matrix(pcol,,1)
    v <- matrix(prow,1) %*% Pbar
  
    delta2 <- 3 * sum(P * P0bar)
    dfdpi <- P0bar + rep(u, ncol) + rep(v, each=nrow)
    sigma <- sqrt(9 * sum(P * dfdpi^2) - 9 * delta2^2)
    dist95 <- 1.959964 * sigma / sqrt(N)
  
    if(abs(delta2) < 1) {
        fdelta2 <- 0.5 * log((1 + delta2) / (1 - delta2))
        fsigma <- sigma/(1-delta2^2)
        fdist95 <- 1.959964*fsigma/sqrt(N)
    }
  
    finv = function(x) 1-2/(exp(2*x) + 1)
    return(c(delta2, delta2 - dist95, delta2 + dist95,
             finv(fdelta2 - fdist95), finv(fdelta2 + fdist95)))    
}
    
condis <- function(P) {
  ## The input P is a matrix of data counts.
  nrow <- dim(P)[1]
  ncol <- dim(P)[2]
  if(nrow<2 | ncol<2) return(0)
  
  N <- sum(P)  # sample size

  ## P is now a probability distribution.
  P <- P/N
  
  ## P0 is the expected distribution under independence.
  prow <- P %*% rep(1,ncol)
  pcol <- rep(1, nrow) %*% P
  P0 <- prow %*% pcol
  
  Pbar <- P0bar <- matrix(,nrow,ncol)
  
  ## For each position (i,j), calculate the probabilities for the 4 sections:
  ## aa1[j] = probability for the section (<i, <j)
  ## aa2[j] = probability for the section (<i, >j)
  ## bb1[j] = probability for the section (>i, <j)
  ## bb2[j] = probability for the section (>i, >j)
  aa <- rep(0,ncol)
  bb <- pcol
  for(i in 1:nrow) {
    if(i > 1)
        aa <- aa + P0[i-1,]
    
    bb <- bb - P0[i,]
    aa1 <- c(0, cumsum(aa[-ncol]))
    aa2 <- sum(aa) - cumsum(aa)
    bb1 <- c(0, cumsum(bb[-ncol]))
    bb2 <- sum(bb) - cumsum(bb)
    P0bar[i,] <- aa1 + bb2 - aa2 - bb1
  }
  
  aa <- rep(0,ncol)
  bb <- pcol
  for(i in 1:nrow) {
    if(i > 1)
        aa <- aa + P[i-1,]

    bb <- bb - P[i,]
    aa1 <- c(0, cumsum(aa[-ncol]))
    aa2 <- sum(aa) - cumsum(aa)
    bb1 <- c(0, cumsum(bb[-ncol]))
    bb2 <- sum(bb) - cumsum(bb)
    Pbar[i,] <- aa1 + bb2 - aa2 - bb1
  }
  
  u <- Pbar %*% matrix(pcol,,1)
  v <- matrix(prow,1) %*% Pbar
  
  delta2 <- 3 * sum(P * P0bar)
  dfdpi <- P0bar + rep(u, ncol) + rep(v, each=nrow)
  sigma <- sqrt(9 * sum(P * dfdpi^2) - 9 * delta2^2)
  dist95 <- 1.959964 * sigma / sqrt(N)
  
  if(abs(delta2) < 1) {
    fdelta2 <- 0.5 * log((1 + delta2) / (1 - delta2))
    fsigma <- sigma/(1-delta2^2)
    fdist95 <- 1.959964*fsigma/sqrt(N)
  }
  
  finv = function(x) 1-2/(exp(2*x) + 1)
  return(c(delta2, delta2 - dist95, delta2 + dist95,
           finv(fdelta2 - fdist95), finv(fdelta2 + fdist95)))
}
