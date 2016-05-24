cvpvs.gaussian <-
function(X, Y, cova = c('standard', 'M', 'sym')) {
  cova <- match.arg(cova)

  # Adjust input
  X <- as.matrix(X)
  n <- NROW(X)
  dimension <- NCOL(X)
  
  Ylevels <- levels(factor(Y))
  Y <- as.integer(factor(Y))
  L <- max(Y)

  # Stop if lengths of X[,1] and Y do not match
  if(length(Y) != length(X[,1])) {
    stop('length(Y) != length(X[,1])')
  }
  
  # Computation of nvec: an L-dimensional column vector, where nvec[b]
  # is the number of training observations belonging to class b.
  nvec <- rep(0, L)
  for(b in seq_len(L)) {
    nvec[b] <- sum(Y == b)
    
    # Stop if there are less than two observations from class b
    if(nvec[b] == 0) {
      stop(paste('no observations from class',
                 as.character(b),'!'))
    }
    if(nvec[b] == 1) {
      stop(paste('only one observation from class',
                 as.character(b),'!'))
    }
  }
  

  # Compute sigma
  switch(cova,
         'standard' = {
           # Compute mu
           mu <- matrix(0, L, dimension)
           for(b in seq_len(L)) {
             mu[b, ] = colMeans(X[Y == b, , drop = FALSE])
           }
           
           sigma <- sigmaSt(X = X, Y = Y, L = L,
                                      dimension = dimension, n = n, mu = mu)
         },
         'M' = {
           tmp <- MVTMLE.LDA(X = X, Y = Y, L = L, n = n, nu = 1,
                                       M0=NULL,B0=NULL,
                                       delta=10^(-7),steps=FALSE)
           mu <- tmp$M
           sigma <- tmp$S
           B <- tmp$B
         },
         'sym' = {
           # Compute mu
           mu <- matrix(0, L, dimension)
           for(b in seq_len(L)) {
             mu[b, ] = colMeans(X[Y == b, , drop = FALSE])
           }
           
           tmp <- sigmaSym(X = X, Y = Y, L = L,
                                     dimension = dimension, n = n,
                                     nvec = nvec,
                                     B = NULL,
                                     nu=0,
                                     delta=10^(-7), prewhitened=TRUE,
                                     steps=FALSE, nmax=500
                                     )
           sigma <- tmp$S
           B <- tmp$B
         } )
  
  sigma.inv <- solve(sigma)
  
  # Computation of the cross-validated P-values:
  PV <- matrix(0, n, L)

  
  # P-values for the correct classes, i.e. PV[i, Y[i]]:
  T <- rep(0, n)

  mu.diff <- NULL
  mu.2 <- NULL
  for(b in seq_len(L - 1)) {
    mu.diff <- cbind(mu.diff, mu[b, ] - t(mu[(b+1):L, , drop = FALSE]))
    mu.2 <- cbind(mu.2, (mu[b, ] + t(mu[(b+1):L, , drop = FALSE])) / 2)
  }
  delta.tmp <- qr.solve(sigma, mu.diff)
  delta.tmpA <- array(0, c(dimension, L, L))
  mu.2A <- array(0, c(dimension, L, L))
  z <- 0
  for(b in seq_len(L - 1)) {
    for(th in (b + 1):L) {
      z <- z + 1
      delta.tmpA[ , b, th] <-  delta.tmp[ , z]
      delta.tmpA[ , th, b] <- -delta.tmp[ , z]
      mu.2A[ , b, th] <- mu.2A[ , th, b] <- mu.2[ , z]
    }
  }
  
  for(i in seq_len(n)) {
    for(b in seq_len(L)[-Y[i]]) {
        T[i] <- T[i] + nvec[b] *
          exp((X[i, ] - mu.2A[ , Y[i], b]) %*% delta.tmpA[ , b, Y[i]]) 
    }
  }
  
  for(th in seq_len(L)) {
    thIndex <- Y == th
    Tv <- -T[thIndex]
    PV[thIndex, th] <- CompPVs(Tv)
  }

    
  # P-values for other classes, i.e. PV[i,th], th != Y[i]:
  for(i in seq_len(n)) {
    for(th in seq_len(L)[-Y[i]]) {
      T <- rep(0, n)

      # Adjust nvec
      nvec.tmp <- nvec
      nvec.tmp[c(Y[i], th)] <- nvec.tmp[c(Y[i], th)] + c(-1, 1)
      
      # Adjust sigma
      switch(cova,
             'standard' = {
               # Adjust mu
               mu.tmp <- mu
               mu.tmp[Y[i], ] <- mu[Y[i], ] - (X[i, ] - mu[Y[i], ]) / (nvec[Y[i]] - 1)
               mu.tmp[th, ] <- mu[th, ] + (X[i, ] - mu[th, ]) / (nvec[th] + 1)
               
               sigma.tmp <- sigma - (tcrossprod(X[i, ] - mu[Y[i], ])
                                       / (1 - 1 / nvec[Y[i]])
                                     - tcrossprod(X[i, ] - mu[th, ])
                                         / (1 + 1 / nvec[th])) / (n - L)
             },
             'M' = {
               Y.tmp <- Y
               Y.tmp[i] <- th
               tmp <- MVTMLE.LDA(X = X, Y = Y.tmp, L = L, n = n, nu = 1,
                                           M0=mu,B0=B,
                                           delta=10^(-7),steps=FALSE)
               mu.tmp <- tmp$M
               sigma.tmp <- tmp$S
             },
             'sym' = {
               # Adjust mu
               mu.tmp <- mu
               mu.tmp[Y[i], ] <- mu[Y[i], ] - (X[i, ] - mu[Y[i], ]) / (nvec[Y[i]] - 1)
               mu.tmp[th, ] <- mu[th, ] + (X[i, ] - mu[th, ]) / (nvec[th] + 1)
               
               Y.tmp <- Y
               Y.tmp[i] <- th
               sigma.tmp <- sigmaSym(X = X, Y = Y.tmp, L = L,
                                               dimension = dimension, n = n,
                                               nvec = nvec.tmp, B = B, nu=0,
                                               delta=10^(-7), prewhitened=TRUE,
                                               steps=FALSE, nmax=500)$S
             } )


      # delta.tmp[b,] := sigma^(-1) (mu[b] - mu[th])
      delta.tmp <- t(qr.solve(sigma.tmp, t(mu.tmp) - mu.tmp[th, ]))
      # mu.th[b,] := (mu.tmp[th, ] + mu.tmp[b, ]) / 2
      mu.th <- t(t(mu.tmp) + mu.tmp[th, ]) / 2
      
      # Compute T
      thIndex <- c(i, which(Y == th))
      for(j in thIndex){
        for(b in seq_len(L)[-th]){
          T[j] <- T[j] + nvec.tmp[b] *
                          exp((X[j, ] - mu.th[b, ]) %*% delta.tmp[b, ])
        }
      }

      Tv <- T[thIndex]
      PV[i, th] <- sum(Tv >= Tv[1]) / nvec.tmp[th]
    }
  }
  dimnames(PV)[[2]] <- Ylevels
  return(PV)
}
