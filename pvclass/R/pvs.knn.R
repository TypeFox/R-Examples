pvs.knn <-
function(NewX, X, Y, k = NULL,
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
  distance <- match.arg(distance)
  
  # Computation of nvec = vector with the numbers of training
  # observations from each class:
  nvec <- rep(0, L)
  for(b in seq_len(L)) {
    nvec[b] <- sum(Y == b)
    
    # Stop if there are less than two observations from class b
    if(nvec[b] == 0) {
      stop(paste('no observations from class', as.character(b),'!'))
    }
    if(nvec[b] == 1) {
      stop(paste('only one observation from class', as.character(b), '!'))
    }
  }

  
  # If k is not a single integer, choose optimal k
  if(length(k) > 1 | is.null(k)) {
    if(is.null(k)) {
      k <- 2:ceiling(length(Y) / 2)
    }
    opt.k <- matrix(0, nr, L)

    switch(distance,
           'euclidean' = {
             # Use Euclidean distance
             di <- as.matrix(dist(rbind(X, NewX)))
             dimnames(di) <- NULL
    
             for(th in seq_len(L)) {
               thIndex <- c(Y == th, TRUE)        
               for(i in seq_len(nr)) {
                 di.tmp <- di[seq_len(n), c(seq_len(n), n+i)]
                 T <- matrix(0, n, length(k))
                 for(j in which(Y != th)) {
                   sdi <- sort(di.tmp[j, ])
                   thIndex[j] <- TRUE
                   for(l in seq_along(k)) {
                     Rl <- sdi[k[l]]
                     T[j, l] <- sum((di.tmp[j, thIndex] <= Rl)) / k[l]
                   }
                   thIndex[j] <- FALSE
                 }
                 T <- colSums(T)
                 opt.k[i, th] <- k[which.min(T)]
               }
             }
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

               di.tmp <- as.matrix(dist(Xtmp))
               
               for(th in seq_len(L)) {
                 T <- matrix(0, n, length(k))
                 thIndex <- c(Y == th, TRUE)
                 for(j in which(Y != th)) {
                   sdi <- sort(di.tmp[j, ])
                   thIndex[j] <- TRUE
                   for(l in seq_along(k)) {
                     Rl <- sdi[k[l]]
                     T[j, l] <- sum((di.tmp[j, thIndex] <= Rl)) / k[l]
                   }
                   thIndex[j] <- FALSE
                 }
                 T <- colSums(T)     
                 opt.k[i, th] <- k[which.min(T)]
               }
             }
           },
           'mahalanobis' = {
             # Use Mahalanobis distance
             di.tmp <- matrix(0, n, n + 1)

             
             # Compute sigma      
             switch(cova,
                    'standard' = {
                      # Compute mu
                      mu <- matrix(0, L, dimension)
                      for(b in seq_len(L)) {
                        mu[b, ] = colMeans(X[Y == b, , drop = FALSE])
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
        
             Nk <- rep(0, n + 1)
        
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
                 
                 T <- matrix(0, n, length(k))
                 for(j in which(Y != th)) {
                   sdi <- sort(di.tmp[j, ])
                   thIndex[j] <- TRUE
                   for(l in seq_along(k)){
                     Rl <- sdi[k[l]]
                     T[j, l] <- sum((di.tmp[j, thIndex] <= Rl)) / k[l]
                   }
                   thIndex[j] <- FALSE
                 }
                 T <- colSums(T)
                 opt.k[i, th] <- k[which.min(T)]
               }
             }
           })

    PV <- matrix(0, nr, L)
    for(j in unique(as.vector(opt.k))) {
      PV[opt.k == j] <- pvs.knn(NewX = NewX, X = X, Y = Y, k = j,
                                distance = distance, cova = cova)[opt.k == j]      
    }
    
    dimnames(opt.k)[[2]] <- dimnames(PV)[[2]] <- Ylevels
    attributes(PV)$opt.k <- opt.k
    dimnames(PV)[[2]] <- Ylevels
    return(PV)
  }
  

  
  
  PV <- matrix(0, nr, L)

  switch(distance,
         'mahalanobis' = {
           # Use Mahalanobis distance
           di <- matrix(0, n + 1, n + 1)
    

           
           # Compute sigma      
           switch(cova,
                  'standard' = {
                    # Compute mu
                    mu <- matrix(0, L, dimension)
                    for(m in seq_len(L)) {
                      mu[m, ] = colMeans(X[Y == m, , drop = FALSE ])
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
           
           Nk <- rep(0, n + 1)
           
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
               for(m in seq_len(n + 1)) {
                 di[m, ] <- mahalanobis(Xtmp, Xtmp[m, ], sigmatmp.inv,
                                        inverted = TRUE)
               }        
        
               sdi <- apply(di, 1, sort)
               Rk <- sdi[k, ]
               Nk <- rowSums(di[thIndex, thIndex] <= Rk[thIndex])
               
               # Compute PV[i, th]
               PV[i, th] <- sum(Nk <= tail(Nk, 1)) / nvec.tmp[th]        
             }
           }
           dimnames(PV)[[2]] <- Ylevels
           return(PV)
         },
         'ddeuclidean' = {  
           # Use data driven Euclidean distance
           Nk <- rep(0, n + 1)
           for(i in seq_len(nr)) {		
             # Add new observation NewX[i, ] to Xtmp
             Xtmp <- rbind(X, NewX[i, ])
             
             # Compute data driven euclidean distances
             for(m in seq_len(dimension)) {
               Xtmp[ , m] <- Xtmp[ , m] / sd(Xtmp[ , m])
             }
             
             di <- as.matrix(dist(Xtmp))
             sdi <- apply(di, 1, sort)
             Rk <- sdi[k, ]
             
             for(th in seq_len(L)) {
               # Add NewX[i] temporarily to group th:
               thIndex <- c(Y == th, TRUE)
               Nk <- rowSums(di[thIndex, thIndex] <= Rk[thIndex])
               
               # Compute PV[i,th]
               PV[i, th] <- sum(Nk <= tail(Nk, 1)) / (nvec[th] + 1)
             }
           }
           dimnames(PV)[[2]] <- Ylevels
           return(PV)
         },
         'euclidean' = {
           # Use Euclidean distance
           di <- as.matrix(dist(rbind(X, NewX)))
           Nk <- matrix(0, n, L + 1)
           Nkm1 <- matrix(0, n, L)
           sdi <- apply(di[seq_len(n),seq_len(n)], 1, sort)
           Rk <- sdi[k, ]
           rkm1 <- sdi[k - 1, ]
           for(j in 1:n) {
             for(th in 1:L) {
               thIndex <- which(Y==th)
               Nk[j,th] <- sum(di[1:n,1:n][j,thIndex]<=Rk[j])
               Nkm1[j,th] <- sum(di[1:n,1:n][j,thIndex]<=rkm1[j])
             }
           }
           
           Nk <- rbind(Nk, c(rep(0, L), 1))
           
           for(i in seq_len(nr)) {		
             Nk_tmp <- Nk
             
             # Add new observation NewX[i, ] to Xtmp
             Xtmp <- rbind(X, NewX[i, ])
             di_tmp <- di[c(seq_len(n), n + i), c(seq_len(n), n + i)]
             sdi_tmp <- apply(di_tmp, 1, sort)
             NewRk <- sdi_tmp[k, ]
             
             # Determine Nk_tmp[n + 1, seq_len(L)]
             for(th in seq_len(L)) {
               thIndex <- which(Y == th)
               Nk_tmp[n + 1, th] <- sum(di_tmp[n + 1, thIndex] <= NewRk[n + 1])
             }
             # Update Nk_tmp[1:n,]
             JJ1 <- which(di_tmp[n + 1, seq_len(n)] < Rk)
             JJ2 <- which(di_tmp[n + 1, seq_len(n)] <= Rk)
             Nk_tmp[JJ1, seq_len(L)] <- Nkm1[JJ1, ]
             Nk_tmp[JJ2, L + 1] <- 1
             
             for(th in seq_len(L)) {
               # Add NewX[i] temporarily to group th:
               thIndex <- c(which(Y == th), n + 1)
               Nk2 <- Nk_tmp[thIndex, th] + Nk_tmp[thIndex, L + 1]
               
               # Compute PV[i, th]
               PV[i, th] <- sum(Nk2 <= tail(Nk2, 1)) / (nvec[th] + 1)
             }
           }
           dimnames(PV)[[2]] <- Ylevels
           return(PV) } )
}
