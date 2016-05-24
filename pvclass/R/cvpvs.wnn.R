cvpvs.wnn <-
function(X, Y, wtype = c('linear', 'exponential'), W = NULL, tau = 0.3,
         distance = c('euclidean', 'ddeuclidean', 'mahalanobis'),
         cova = c('standard', 'M', 'sym')) {

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

  
  wtype <- match.arg(wtype)
  cova <- match.arg(cova)
  distance <- match.arg(distance)

  # Computation of nvec: an L-dimensional column vector, where nvec[b]
  # is the number of training observations belonging to class b.
  nvec <- rep(0, L)
  for(b in seq_len(L)) {
    nvec[b] <- sum(Y == b)
    
    # Stop if there are less than two observations from class b
    if(nvec[b] == 0) {
      stop(paste('no observations from class',
                 as.character(b), '!'))
    }
    if(nvec[b] == 1) {
      stop(paste('only one observation from class',
                 as.character(b), '!'))
    }
  }
  
  # If tau is not a single value, search for the optimal tau.
  if( is.null(W) & (length(tau) > 1 | is.null(tau)[1]) ) {
    switch(wtype,
           # Use exponential weights
           'exponential' = {
             if(is.null(tau)[1]) { tau <- c(1, 5, 10, 20) }
             W <- matrix(1 - seq_len(n) / n, length(tau),
                         n, byrow = TRUE)^tau
           },
           # Use linear weights
           'linear' = {
             if(is.null(tau)[1]) { tau <- seq(0.1, 0.9, 0.1) }
             W <- pmax(1 - matrix(seq_len(n) / n, length(tau), n, byrow = TRUE)
                       / tau, 0)
           } )
    
    Wsum <- rowSums(W)
    
    opt.tau <- matrix(0, n, L)

    # optimal tau for correct classes
    switch(distance,
           'euclidean' = {
             # Use Euclidean distance
             di.tmp <- as.matrix(dist(X))
             
             for(th in seq_len(L)) {
               T <- matrix(0, n, length(tau))
               WXtmp <- matrix(0, length(tau), n)
               thIndex <- Y == th
               for(j in which(Y != th)) {
                 thIndex[j] <- TRUE
                 for(l in seq_along(tau)) {
                   WXtmp[l, order(di.tmp[j, ])] <- W[l, ]
                 }
                 T[j, ] <- (WXtmp %*% thIndex) / Wsum
                 thIndex[j] <- FALSE
               }
               T <- colSums(T)
               opt.tau[Y == th, th] <- tau[which.min(T)]
             }
           },
           'ddeuclidean' = {
             # Use data driven Euclidean distance
             for(m in seq_len(dimension)) {
               X[ , m] <- X[ , m] / sd(X[ , m])
             }
             di.tmp <- as.matrix(dist(X))
             
             for(th in seq_len(L)) {
               T <- matrix(0, n, length(tau))
               WXtmp <- matrix(0, length(tau), n)
               thIndex <- Y == th
               for(j in which(Y != th)) {
                 thIndex[j] <- TRUE
                 for(l in seq_along(tau)) {
                   WXtmp[l, order(di.tmp[j, ])] <- W[l, ]
                 }
                 T[j, ] <- (WXtmp %*% thIndex) / Wsum
                 thIndex[j] <- FALSE
               }
               T <- colSums(T)
               opt.tau[Y == th, th] <- tau[which.min(T)]
             }
           },
           'mahalanobis' = {
             # Use Mahalanobis distance
             di.tmp <- matrix(0, n, n)
    
             # Compute mu
             mu <- matrix(0, L, dimension)
             for(m in seq_len(L)) {
               mu[m, ] = colMeans(X[Y == m, , drop = FALSE])
             }
             
             # Compute sigma      
             switch(cova,
                    'standard' = {
                      sigma <- sigmaSt(X = X, Y = Y, L = L,
                                                 dimension = dimension, n = n,
                                                 mu = mu)
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

             # Compute Mahalanobis distances
             sigma.inv <- solve(sigma)
             for(m in seq_len(n)) {
               di.tmp[m, ] <- mahalanobis(X, X[m, ], sigma.inv, inverted = TRUE)
             }
             
             for(th in seq_len(L)) {
               T <- matrix(0, n, length(tau))
               WXtmp <- matrix(0, length(tau), n)
               thIndex <- Y == th
               for(j in which(Y != th)) {
                 thIndex[j] <- TRUE
                 for(l in seq_along(tau)) {
                   WXtmp[l, order(di.tmp[j, ])] <- W[l, ]
                 }
                 T[j, ] <- (WXtmp %*% thIndex) / Wsum
                 thIndex[j] <- FALSE
               }
               T <- colSums(T)
               opt.tau[Y == th, th] <- tau[which.min(T)]          
             }
           } )
    
    
    # optimal tau other classes
    if(distance != 'mahalanobis') {
      # Use Euclidean distance or data driven Euclidean distance
      for(th in seq_len(L)) {
        thIndex <- Y == th
        for(i in which(Y != th)) {
          T <- matrix(0, n, length(tau))
          thIndex[i] <- TRUE
          for(j in which(!thIndex)) {
            thIndex.tmp <- thIndex
            thIndex.tmp[j] <- TRUE
            for(l in seq_along(tau)) {
              WXtmp[l, order(di.tmp[j, ])] <- W[l, ]
            }
            T[j, ] <- (WXtmp %*% thIndex) / Wsum
          }
          
          T <- colSums(T)
          opt.tau[i, th] <- tau[which.min(T)]
          thIndex[i] <- FALSE
        }
      }
    } else {
      # Use Mahalanobis distance
      for(th in seq_len(L)) {
        thIndex <- Y == th
        for(i in which(Y != th)) {
          # Add NewX[i] temporarily to group th:
          Y.tmp <- Y
          Y.tmp[i] <- th
          thIndex <- Y.tmp == th
          nvec.tmp <- nvec
          nvec.tmp[c(Y[i], th)] <- c(nvec[Y[i]] - 1, nvec[th] + 1)
         
          # Adjust sigma
          switch(cova,
                 'standard' = {
                   sigma.tmp <- sigma - (tcrossprod(X[i,]-mu[Y[i],]) /
                                         (1 - 1 / nvec[Y[i]]) -
                                         tcrossprod(X[i,]-mu[th,]) /
                                         (1 - 1 / nvec[th]) ) / (n-L)
                 },
                 'M' = {
                   sigma.tmp <- MVTMLE.LDA(X = X, Y = Y.tmp, L = L,
                                                     n = n, nu = 1,
                                                     M0=mu,B0=B,
                                                     delta=10^(-7),steps=FALSE)$S
                 },
                 'sym' = {
                   sigma.tmp <- sigmaSym(X = X, Y = Y.tmp, L = L,
                                                   dimension = dimension, n = n,
                                                   nvec = nvec.tmp, B = B, nu=0,
                                                   delta=10^(-7), prewhitened=TRUE,
                                                   steps=FALSE, nmax=500)$S
                 } )
          
          # Compute Mahalanobis distances
          sigmatmp.inv <- solve(sigma.tmp)
          for(m in seq_len(n)) {
            di.tmp[m, ] <- mahalanobis(X, X[m, ], sigmatmp.inv, inverted = TRUE)
          }
          
          T <- matrix(0, n, length(tau))
          thIndex[i] <- TRUE
          for(j in which(!thIndex)) {
            thIndex[j] <- TRUE
            for(l in seq_along(tau)) {
              WXtmp[l, order(di.tmp[j, ])] <- W[l, ]
            }
            T[j, ] <- (WXtmp %*% thIndex) / Wsum
            thIndex[j] <- FALSE
          }
          
          T <- colSums(T)
          opt.tau[i, th] <- tau[which.min(T)]          
          thIndex[i] <- FALSE
        }
      }
    }
    
    PV <- matrix(0, n, L)
    for(j in unique(c(opt.tau))) {
      PV[opt.tau == j] <- cvpvs.wnn(X = X, Y = Y, W = W[which(tau == j)[1], ],
       distance = distance, cova = cova)[opt.tau == j]
    }
    
    attributes(PV)$opt.tau <- opt.tau
    dimnames(PV)[[2]] <- Ylevels
    return(PV)
  }
  
 
  # Computation of weight function W
  if(!is.null(W)) {
    if(length(W) != n) {
      stop(paste('length(W) != length(Y)'))
    }
  } else {
    switch(wtype,
           # Use exponential weights
           'exponential' = W <- matrix(1 - seq_len(n) / n, length(tau), n,
                                 byrow = TRUE)^tau,
           # Use linear weights
           'linear' = W <- pmax(1 - seq_len(n) / n / tau, 0))
  }

    
  WX <- matrix(0, n, n)
  
  if(distance == "mahalanobis") {
    # Use Mahalanobis distance
    
    # Compute mu
    mu <- matrix(0, L, dimension)
    for(m in seq_len(L)) {
      mu[m, ] = colMeans(X[Y == m, , drop = FALSE])
    }
    
    # Compute sigma      
    switch(cova,
           'standard' = {
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



    # Compute Mahalanobis distances
    sigma.inv <- solve(sigma)
    for(i in seq_len(n)) {
      WX[i, order(mahalanobis(X, X[i, ], sigma.inv, inverted = TRUE))] <- W
    }
        
    # Computation of the cross-validated P-values:
    PV <- matrix(0, n, L)

    # P-values for the correct classes, i.e. PV[i, Y[i]]:
    for(th in seq_len(L)) {
      thIndex <- Y == th
      Tv <- WX[thIndex, ] %*% thIndex
      PV[thIndex, th] <- CompPVs(Tv)
    }

    # P-values for other classes, i.e. PV[i, th], th != Y[i]:
    for(i in seq_len(n)) {
      for(th in seq_len(L)[-Y[i]]) {
        # Adjust nvec
        nvec.tmp <- nvec
        nvec.tmp[c(Y[i], th)] <- nvec.tmp[c(Y[i], th)] + c(- 1, 1)
        
        # Adjust sigma
        switch(cova,
               'standard' = {
                 sigma.tmp <- sigma - (tcrossprod(X[i,]-mu[Y[i],]) /
                                         (1 - 1 / nvec[Y[i]]) -
                                       tcrossprod(X[i,]-mu[th,]) /
                                           (1 - 1 / nvec[th]) ) / (n-L)
               },
               'M' = {
                 Y.tmp <- Y
                 Y.tmp[i] <- th
                 sigma.tmp <- MVTMLE.LDA(X = X, Y = Y.tmp, L = L, n = n, nu = 1,
                                                   M0=mu,B0=B,
                                                   delta=10^(-7),steps=FALSE)$S
               },
               'sym' = {
                 Y.tmp <- Y
                 Y.tmp[i] <- th
                 sigma.tmp <- sigmaSym(X = X, Y = Y.tmp, L = L,
                                                 dimension = dimension, n = n,
                                                 nvec = nvec.tmp, B = B, nu=0,
                                                 delta=10^(-7), prewhitened=TRUE,
                                                 steps=FALSE, nmax=500)$S
               } )
        
        thIndex <- c(i, which(Y == th))
        WX <- matrix(0, n, n)

        # Compute Mahalanobis distances
        sigmatmp.inv <- solve(sigma.tmp)
        for(j in thIndex) {
          WX[j, order(mahalanobis(X, X[j, ], sigmatmp.inv, inverted = TRUE))] <- W
        }
        
        HH <- rep(0, n)
        HH[thIndex] <- 1

        Tv <- WX[thIndex, ] %*% HH
        PV[i, th] <- sum(Tv <= Tv[1]) / nvec.tmp[th]
      }
    }   
    
    dimnames(PV)[[2]] <- Ylevels
    return(PV)
    
  } else {
    # Use Euclidean distance           
    if(distance == 'ddeuclidean') {
      # Use data driven Euclidean distance
      for(i in seq_len(dimension)) {
        X[ , i] <- X[ , i] / sd(X[ , i])
      }
    }
    
    Distance <- as.matrix(dist(X))
    for(i in seq_len(n)) {
      WX[i, order(Distance[i, ])] <- W
    }
    
    # Computation of the cross-validated P-values:
    PV <- matrix(0, n, L)
  
    # P-values for the correct classes, i.e. PV(i, Y[i]):
    for(th in seq_len(L)) {
      thIndex <- Y == th
      Tv <- WX[thIndex, ] %*% thIndex
      
      PV[thIndex, th] <- CompPVs(Tv)
    }
  
    # P-values for other classes, i.e. PV(i, th), th != Y(i):
    for(i in seq_len(n)) {
      for(th in seq_len(L)[-Y[i]]) {
        # Add observation i temporarily to class th
        thIndex <- c(i, which(Y == th))
        HH <- rep(0, n)
        HH[thIndex] <- 1
        
        Tv <- WX[thIndex, ] %*% HH 
        PV[i, th] <- sum(Tv <= Tv[1]) / (nvec[th] + 1)
      }		
    }
    dimnames(PV)[[2]] <- Ylevels
    return(PV)
  }
}
