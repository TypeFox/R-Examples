cvpvs.knn <-
function(X, Y, k = NULL,
         distance = c('euclidean', 'ddeuclidean', 'mahalanobis'),
         cova = c('standard', 'M', 'sym')) {
  
  cova <- match.arg(cova)
  distance <- match.arg(distance)
  
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

  if(any(k >= n)) stop('k >= sample size!')
  
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
  
  # If k is not a single integer, search for the optimal k.
  if(length(k) > 1 | is.null(k)) {
    if(is.null(k)) {
      k <- 2:ceiling(length(Y) / 2)
    }
    opt.k <- matrix(0, n, L)
    
    # optimal k for correct classes
    switch(distance,
           'euclidean' = {
             # Use Euclidean distance
             di.tmp <- as.matrix(dist(X))
             
             for(th in seq_len(L)) {
               T <- matrix(0, n, length(k))
               thIndex <- Y == th
               for(j in which(Y != th)) {
                 sdi <- sort(di.tmp[j, ])
                 thIndex[j] <- TRUE
                 for(l in seq_along(k)) {
                   Rl <- sdi[k[l]]
                   T[j, l] <- sum(di.tmp[j, thIndex] <= Rl) / k[l]
                 }
               }
               thIndex[j] <- FALSE

               T <- colSums(T)
               opt.k[thIndex, th] <- k[which.min(T)]
             }
           },
           'ddeuclidean' = {
             # Use data driven Euclidean distance
             for(m in seq_len(dimension)) {
               X[ , m] <- X[ , m] / sd(X[ , m])
             }
             di.tmp <- as.matrix(dist(X))
             
             for(th in seq_len(L)) {
               T <- matrix(0, n, length(k))
               thIndex <- Y == th
               for(j in which(Y != th)) {
                 sdi <- sort(di.tmp[j, ])
                 thIndex[j] <- TRUE
                 for(l in seq_along(k)) {
                   Rl <- sdi[k[l]]
                   T[j, l] <- sum(di.tmp[j, thIndex] <= Rl) / k[l]
                 }
               }
               thIndex[j] <- FALSE
               
               T <- colSums(T)
               opt.k[thIndex, th] <- k[which.min(T)]          
             }
           },
           'mahalanobis' = {
             # Use Mahalanobis distance
             di.tmp <- matrix(0, n, n)
    

             
             # Compute sigma      
             switch(cova,
                    'standard' = {
                      # Compute mu
                      mu <- matrix(0, L, dimension)
                      for(m in seq_len(L)) {
                        mu[m, ] = colMeans(X[Y == m, , drop = FALSE])
                      }
                      sigma <- sigmaSt(X = X, Y = Y, L = L,
                                                 dimension = dimension,
                                                 n = n, mu = mu)
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
               T <- matrix(0, n, length(k))
               thIndex <- Y == th
               for(j in which(Y != th)) {
                 sdi <- sort(di.tmp[j, ])
                 thIndex[j] <- TRUE
                 for(l in seq_along(k)) {
                   Rl <- sdi[k[l]]
                   T[j, l] <- sum(di.tmp[j, thIndex] <= Rl) / k[l]
                 }
               }
               thIndex[j] <- FALSE

               T <- colSums(T)
               opt.k[thIndex, th] <- k[which.min(T)]                            
             }

           } )
         

    
    
    # optimal k other classes
    if(distance == "euclidean" | distance == 'ddeuclidean') {
      # Use Euclidean distance or data driven Euclidean distance
      for(th in seq_len(L)) {
        thIndex <- Y == th
        for(i in which(Y != th)) {
          T <- matrix(0, n, length(k))
          thIndex[i] <- TRUE
          for(j in which(!thIndex)) {
            sdi <- sort(di.tmp[j, ])
            thIndex[j] <- TRUE
            for(l in seq_along(k)) {
              Rl <- sdi[k[l]]
              T[j, l] <- sum(di.tmp[j, thIndex] <= Rl) / k[l]
            }
            thIndex[j] <- FALSE
          }
          thIndex[i] <- FALSE
          
          T <- colSums(T)
          opt.k[i, th] <- k[which.min(T)]
        }
      }
    } else {
      # Use Mahalanobis distance
      for(th in seq_len(L)) {
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
                   sigma.tmp <- MVTMLE.LDA(X = X, Y = Y, L = L, n = n, nu = 1,
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
            di.tmp[m, ] <- mahalanobis(X, X[m, ], sigmatmp.inv,
                                       inverted = TRUE)
          }

          T <- matrix(0, n, length(k))
          thIndex.tmp <- thIndex
          for(j in which(!thIndex)) {
            sdi <- sort(di.tmp[j, ])
            thIndex.tmp[j] <- TRUE
            for(l in seq_along(k)) {
              Rl <- sdi[k[l]]
              T[j, l] <- sum(di.tmp[j, thIndex.tmp] <= Rl) / k[l]
            }
          }
          thIndex.tmp[j] <- FALSE
          
          T <- colSums(T)
          opt.k[i, th] <- k[which.min(T)]
        }
      }
    }
    PV <- matrix(0, n, L)

    for(j in unique(c(opt.k))) {
      PV[opt.k == j] <- cvpvs.knn(X = X, Y = Y, k = j, distance = distance,
                                  cova = cova)[opt.k == j]
    }
        
    attributes(PV)$opt.k <- opt.k
    dimnames(PV)[[2]] <- Ylevels
    return(PV)
  }
  
   
  

  
  # Compute distances
  if(distance == "mahalanobis") {
    # Use Mahalanobis distance
   
    # Compute sigma      
    switch(cova,
           'standard' = {
             # Compute mu
             mu <- matrix(0, L, dimension)
             for(m in seq_len(L)) {
               mu[m, ] = colMeans(X[Y == m, , drop = FALSE])
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
    di <- matrix(0, n, n)
    for(m in seq_len(n)) {
      di[m, ] <- mahalanobis(X, X[m, ], sigma.inv, inverted = TRUE)
    }
    
    # Computation of Rk :  a column vector, where Rk[i] is the
    #                      distance between X[i,] and its k-th nearest
    #                      neighbor.
    #                Nk :  Nk[i,b] is the number of training
    #                      observations from class b among the k
    #                      nearest neighbors of X[i,]
    
    Rk <- rep(0, n)
    Nk <- matrix(0, n, L)
    
    sdi <- apply(di, 1, sort)
    Rk <- sdi[k, ]
   
    for(b in seq_len(L)) {
      bIndex <- Y == b
      # get number of observations for each class among k
      # neighbours:
      Nk[ , b] <- rowSums(di[ , bIndex] <= Rk)
    }

    
    # Computation of the cross-validated P-values:
    PV <- matrix(0, n, L)
    
    # P-values for the correct classes, i.e. PV[i, Y[i]]:
    for(th in seq_len(L)) {
      thIndex <- Y == th
      Tv <- Nk[thIndex, th]
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
                 'standard' = sigma.tmp <- sigma -
                   (tcrossprod(X[i, ] - mu[Y[i], ]) / (1 - 1 / nvec[Y[i]])
                      - tcrossprod(X[i, ] - mu[th, ]) / (1 - 1 / nvec[th])
                    ) / (n - L),
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

          # Compute Mahalanobis distances
          sigmatmp.inv <- solve(sigma.tmp)
          for(m in thIndex) {
            di[m, ] <- mahalanobis(X, X[m, ], sigmatmp.inv, inverted = TRUE)
          }

          # Computation of Rk :  a column vector, where Rk[i] is the
          #                      distance between X[i,] and its k-th nearest
          #                      neighbor.
          #                Nk :  Nk[i,b] is the number of training
          #                      observations from class b among the k
          #                      nearest neighbors of X[i,]
          Rk <- rep(0, n)
          Nk <- matrix(0, n, L)
          
          sdi <- apply(di, 1, sort)
          Rk <- sdi[k, ]
          Nk[ , th] <- rowSums(di[ , thIndex] <= Rk)
          
          Tv <- Nk[thIndex, th]
          PV[i, th] <- sum(Tv <= Tv[1]) / nvec.tmp[th]        
      }
    }
    dimnames(PV)[[2]] <- Ylevels
    return(PV)
  } else {
    # Use Euclidean distance

    if(distance == "ddeuclidean") {
      # Use data driven Euclidean distance
      for(i in seq_len(dimension)) {
        X[ , i] <- X[ , i] / sd(X[ , i])
      }
    }
    di  <- as.matrix(dist(X))
  }
  
  # Computation of Rk :  a column vector, where Rk[i] is the
  #                      distance between X[i,] and its k-th nearest
  #                      neighbor.
  #                Nk :  Nk[i,b] is the number of training
  #                      observations from class b among the k
  #                      nearest neighbors of X[i,]
  Rk <- rep(0, n)
  Nk <- matrix(0, n, L)
  
  sdi <- apply(di, 1, sort)
  Rk <- sdi[k, ]
  for(b in seq_len(L)) {
    bIndex <- Y == b
    # get number of observations for each class among k
    # neighbours:
    Nk[ , b] <- rowSums(di[ , bIndex] <= Rk)
  }
         
  # Computation of the cross-validated P-values:
  PV <- matrix(0, n, L)
    
  # P-values for the correct classes, i.e. PV[i, Y[i]]:
  for(th in seq_len(L)) {
    thIndex <- Y == th
    Tv <- Nk[thIndex, th]
    PV[thIndex, th] <- CompPVs(Tv)
  }
    
    
  # P-values for other classes, i.e. PV[i, th], th != Y[i]:
  for(i in seq_len(n)) {
    for(th in seq_len(L)[-Y[i]]) {
      thIndex <- c(i, which(Y == th))
      Nk_tmp  <- Nk
      #nvec_tmp <- nvec
      RkInd <- di[thIndex, i] <= Rk[thIndex]
      Nk_tmp[thIndex, th] <- Nk_tmp[thIndex, th] + RkInd
      
      Tv <- Nk_tmp[thIndex, th]
      PV[i, th] <- sum(Tv <= Tv[1]) / (nvec[th] + 1 )
    }
  }

  dimnames(PV)[[2]] <- Ylevels
  return(PV)
}
