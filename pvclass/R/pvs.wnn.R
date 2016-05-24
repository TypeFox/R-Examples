pvs.wnn <-
function(NewX, X, Y, wtype = c('linear', 'exponential'), W = NULL, tau = 0.3,
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

  NewX <- as.matrix(NewX)
  
  if(dimension > 1 & NCOL(NewX) == 1) {
    NewX <- t(NewX)
  }
  
  if(dimension == 1 & NCOL(NewX) > 1) {
    NewX <- t(NewX)
  }
  
  nr <- NROW(NewX)
  s <- NCOL(NewX)

  # Stop if dimensions of NewX[i,] and X[j,] do not match
  if(s != dimension) {
    stop('dimensions of NewX[i,] and X[j,] do not match!')
  }
  
  cova <- match.arg(cova)
  wtype <- match.arg(wtype)
  distance <- match.arg(distance)

  # If tau is not a single value choose optimal tau
  if( is.null(W) & (length(tau) > 1 | is.null(tau)) ) {
    switch(wtype,
           'exponential' = {
             # Use exponential weights
             if(is.null(tau)) {
               tau <- c(1, 5, 10, 20)
             }
             W <- matrix(1 - seq_len(n + 1) / (n + 1), length(tau), n + 1,
                         byrow = TRUE)^tau
           },
           'linear' = {
             # Use linear weights
             if(is.null(tau)) {
               tau <- seq(0.1, 0.9, 0.1)
             }
             W <- pmax(1 - matrix(seq_len(n + 1) / (n + 1), length(tau), n + 1,
                                  byrow = TRUE) / tau, 0)
           } )

    Wsum <- rowSums(W)

    opt.tau <- matrix(0, nr, L)
    
    switch(distance,
           'euclidean' = {
             # Use Euclidean distance
             di <- as.matrix(dist(rbind(X, NewX)))
             for(i in seq_len(nr)) {
               di.tmp <- di[seq_len(n), c(seq_len(n), n + i)]
               for(th in seq_len(L)) {
                 T <- matrix(0, n, length(tau))
                 WXtmp <- matrix(0, length(tau), n + 1)
                 thIndex <- c(Y == th, TRUE)                 
                 for(j in which(Y != th)) {
                   thIndex[j] <- TRUE
                   for(l in seq_along(tau)) {
                     WXtmp[l, order(di.tmp[j, ])] <- W[l, ]
                   }
                   T[j, ] <- (WXtmp %*% thIndex) / Wsum
                   thIndex[j] <- FALSE
                 }
                 T <- colSums(T)
                 opt.tau[i, th] <- tau[which.min(T)]
               }
             }
           },
           'ddeuclidean' = {
             # Use data driven Euclidean distance
             for(i in seq_len(nr)) {
               # Add new observation NewX[i, ] to Xtmp
               Xtmp <- rbind(X,NewX[i, ])
          
               # Compute data driven euclidean distances
               for(m in seq_len(dimension)) {
                 Xtmp[ , m] <- Xtmp[ , m] / sd(Xtmp[ , m])
               }
               di.tmp <- as.matrix(dist(Xtmp))
               
               for(th in seq_len(L)) {
                 T <- matrix(0, n, length(tau))
                 WXtmp <- matrix(0, length(tau), n + 1)
                 thIndex <- c(Y==th, TRUE)
                 for(j in which(Y != th)) {
                   thIndex[j] <- TRUE
                   for(l in seq_along(tau)) {
                     WXtmp[l, order(di.tmp[j, ])] <- W[l, ]
                   }
                   T[j, ] <- (WXtmp %*% thIndex) / Wsum
                   thIndex[j] <- FALSE
                 }
                 T <- colSums(T)
                 opt.tau[i, th] <- tau[which.min(T)]
               }
             }
           },
           'mahalanobis' = {
             # Use Mahalanobis distance
             di.tmp <- matrix(0, n, n + 1)

             # Computation of nvec = vector with the numbers of training
             # observations from each class:
             nvec <- rep(0, L)
             for(b in seq_len(L)) {
               nvec[b] <- sum(Y == b)
               
               #  Stop if there are less than two observations from class b
               if(nvec[b] == 0) {
                 stop(paste('no observations from class', as.character(b), '!'))
               }
               if(nvec[b]==1) {
                 stop(paste('only one observation from class', as.character(b),
                            '!'))
               }
             }
             

    
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
             
             
             for(i in seq_len(nr)) {		
               # Add new observation NewX[i, ] to Xtmp
               Xtmp <- rbind(X, NewX[i, ])
               
               for(th in seq_len(L)) {
                 # Add NewX[i] temporarily to group th:
                 Ytmp <- c(Y, th)
                 thIndex <- Ytmp == th
                 nvec.tmp <- nvec
                 nvec.tmp[th] <- nvec.tmp[th] + 1
                 
                 # Adjust sigma
                 switch(cova,
                        'standard' = {
                           sigmatmp <- ((n - L) * sigma + 1 / (1 + 1 / nvec[th]) *
                                     tcrossprod(NewX[i, ] - mu[th])
                                     ) / (n + 1 - L)
                         },
                        'M' = {
                          sigmatmp <- MVTMLE.LDA(X = Xtmp, Y = Ytmp, L = L, n = n + 1, nu = 1,
                                                           M0=mu,B0=B,
                                                           delta=10^(-7),steps=FALSE)$S
                        },
                        'sym' = {
                          sigmatmp <- sigmaSym(X = Xtmp, Y = Ytmp, L = L,
                                                          dimension = dimension, n = n + 1,
                                                          nvec = nvec.tmp, B = B, nu=0,
                                                          delta=10^(-7), prewhitened=TRUE,
                                                          steps=FALSE, nmax=500)$S
                        } )
                 
                 # Compute Mahalanobis distances
                 sigmatmp.inv <- solve(sigmatmp)
                 for(m in seq_len(n)) {
                   di.tmp[m, ] <- mahalanobis(Xtmp, Xtmp[m, ], sigmatmp.inv,
                                              inverted = TRUE)
                 }
                 
                 T <- matrix(0, n, length(tau))
                 WXtmp <- matrix(0, length(tau), n + 1)
                 for(j in which(Y != th)) {
                   thIndex[j] <- TRUE
                   for(l in seq_along(tau)) {
                     WXtmp[l, order(di.tmp[j, ])] <- W[l, ]
                   }
                   
                   T[j, ] <- (WXtmp %*% thIndex) / Wsum
                   thIndex[j] <- FALSE
                 }
                 T <- colSums(T)
                 opt.tau[i, th] <- tau[which.min(T)]
               }
             }
           } )
    
    PV <- matrix(0, nr, L)
    for(i in seq_len(nr)) {
      for(j in unique(as.vector(opt.tau))) {
        PV[opt.tau == j] <- pvs.wnn(NewX = NewX, X = X, Y = Y,
                       W = W[which(tau == j)[1], ], distance = distance,
                       cova = cova)[opt.tau == j]
      }
    }
    
    dimnames(opt.tau)[[2]] <- Ylevels
    dimnames(PV)[[2]] <- Ylevels
    attributes(PV)$opt.tau <- opt.tau
    return(PV)
  }
  
  
  
  
  PV <- matrix(0, nr, L)
  
  
  # Computation of nvec = vector with the numbers of training
  # observations from each class:
  nvec <- rep(0, L)
  for(b in seq_len(L)) {
    nvec[b] <- sum(Y == b)
    
    # Stop if there are less than two observations from class b
    if(nvec[b] == 0) {
      stop(paste('no observations from class', as.character(b), '!'))
    }
    if(nvec[b] == 1) {
      stop(paste('only one observation from class', as.character(b),
                 '!'))
    }
  }
  
  # Computation of the weight function W
  if(!is.null(W)) {
    if(length(W) != (n + 1)) {
      stop(paste('length(W) != (length(Y) + 1)'))
    }
  } else {
    switch(wtype,
           # Use exponential weights
           'exponential' = W <- (1 - seq_len(n + 1) / (n + 1))^tau,
           # Use linear weights
           'linear' = W <- pmax(1 - seq_len(n + 1) / (n + 1) / tau, 0))
  }
  
  switch(distance,
         'mahalanobis' = {
           # Use Mahalanobis distance
           WXtmp <- matrix(0, n + 1, n + 1)
           di <- matrix(0, n + 1, n + 1)
           

           
           # Compute sigma
           switch(cova,
                  'standard' = {
                    # Compute mu
                    mu <- matrix(0, L, dimension)
                    for(m in seq_len(L)) {
                      mu[m, ] = colMeans(X[Y == m, , drop = FALSE])
                    }
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
           
           for(i in seq_len(nr)) {		
             # Add new observation NewX[i, ] to Xtmp
             Xtmp <- rbind(X, NewX[i, ])
      
             for(th in seq_len(L)) {
               # Add NewX[i] temporarily to group th:
               Ytmp <- c(Y, th)
               thIndex <- Ytmp == th

               nvec.tmp <- nvec
               nvec.tmp[th] <- nvec.tmp[th] + 1
        
               # Adjust sigma
               switch(cova,
                      'standard' = {
                        sigmatmp <- ((n - L) * sigma + 1 / (1 + 1 / nvec[th]) *
                                     tcrossprod(NewX[i, ] - mu[th])
                                     ) / (n + 1 - L)
                      },
                      'M' = {
                        sigmatmp <- MVTMLE.LDA(X = Xtmp, Y = Ytmp, L = L,
                                                         n = n + 1, nu = 1,
                                                         M0=mu,B0=B,
                                                         delta=10^(-7),steps=FALSE)$S
                      },
                      'sym' = {
                        sigmatmp <- sigmaSym(X = Xtmp, Y = Ytmp, L = L,
                                                       dimension = dimension, n = n + 1,
                                                       nvec = nvec.tmp, B = B, nu=0,
                                                       delta=10^(-7), prewhitened=TRUE,
                                                       steps=FALSE, nmax=500)$S
                      } )
               
               # Compute Mahalanobis distances
               sigmatmp.inv <- solve(sigmatmp)
               for(m in seq_len(n + 1)) {
                 di[m, ] <- mahalanobis(Xtmp, Xtmp[m, ], sigmatmp.inv,
                                            inverted = TRUE)
               }
               
               # Assign weights
               for(j in seq_len(n+1)) {
                 WXtmp[j, order(di[j, ])] <- W
               }
               
               # Compute PV[i, th]
               Tv <- WXtmp[thIndex, ] %*% thIndex
               PV[i, th] <- sum(Tv <= tail(drop(Tv), 1)) / nvec.tmp[th]
             }
           }
           dimnames(PV)[[2]] <- Ylevels
           return(PV)
         },
         'ddeuclidean' = {
           # Use data driven Euclidean distance
           for(i in seq_len(nr)) {		
             # Add new observation NewX[i, ] to Xtmp
             Xtmp <- rbind(X, NewX[i, ])
        
             # Compute data driven euclidean distances
             for(m in seq_len(dimension)) {
               Xtmp[ , m] <- Xtmp[ , m] / sd(Xtmp[ , m])
             }
             di <- as.matrix(dist(Xtmp))
             
             for(th in seq_len(L)) {
               # Add NewX[i] temporarily to group th:
               thIndex <- c(Y == th, TRUE)
               
               # Assign weights
               WXtmp <- matrix(0, n + 1, n + 1)
               for(j in which(thIndex)) {
                 WXtmp[j, order(di[j, ])] <- W
               }
               
               # Compute PV[i, th]
               Tv <- WXtmp[thIndex, ] %*% thIndex
               PV[i, th] <- sum(Tv <= tail(drop(Tv),1)) / (nvec[th] + 1)
             }
           }
           dimnames(PV)[[2]] <- Ylevels
           return(PV)
           
         },
         'euclidean' = {
           # Use Euclidean distance
           di <- as.matrix(dist(rbind(X, NewX)))
           for(i in seq_len(nr)) {		
             # Add new observation NewX[i, ] to Xtmp
             Xtmp <- rbind(X, NewX[i, ])
             di.tmp <- di[c(seq_len(n), n + i), c(seq_len(n), n + i)]
             for(th in seq_len(L)) {
               # Add NewX[i] temporarily to group th:
               thIndex <- c(Y == th, TRUE)
               
               # Assign weights
               WXtmp <- matrix(0, n + 1, n + 1)
               for(j in which(thIndex)) {
                 WXtmp[j, order(di.tmp[j, ])] <- W
               }
               
               # Compute PV[i, th]
               Tv <- WXtmp[thIndex, ] %*% thIndex
               PV[i, th] <- sum(Tv <= tail(drop(Tv), 1)) / (nvec[th] + 1)
             }
           }
           dimnames(PV)[[2]] <- Ylevels
           return(PV)
         } )
}
