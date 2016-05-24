pvs.gaussian <-
function(NewX, X, Y, cova = c('standard', 'M', 'sym')){
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
      stop(paste('only one observation from class' , as.character(b),
                 '!'))
    }
  }
  

  cova <- match.arg(cova)


  
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
  

  for(i in seq_len(nr)) {		
    #Add new observation NewX[i, ] to Xtmp
    X.tmp <- rbind(X, NewX[i, ])

    for(th in seq_len(L)) {
      # Add NewX[i] temporarily to group th:
      Y.tmp <- c(Y, th)
      thIndex <- which(Y.tmp == th)

      nvec.tmp <- nvec
      nvec.tmp[th] <- nvec.tmp[th] + 1

      
      #Adjust sigma
      switch(cova,
             'standard' = {
               # Adjust mu
               mu.tmp <- mu
               mu.tmp[th, ] <- mu[th, ] + (X.tmp[n+1, ] - mu[th, ]) / nvec.tmp[th]

               sigma.tmp <- ((n-L)*sigma + tcrossprod(NewX[i,]-mu[th])
                             / (1 + 1 / nvec[th]) ) / (n + 1 - L)
             },
             'M' = {
               tmp <- MVTMLE.LDA(X = X.tmp, Y = Y.tmp, L = L, n = n + 1, nu = 1,
                                           M0=mu,B0=B,
                                           delta=10^(-7),steps=FALSE)
               mu.tmp <- tmp$M
               sigma.tmp <- tmp$S                      
             },
             'sym' = {
               # Adjust mu
               mu.tmp <- mu
               mu.tmp[th, ] <- mu[th, ] + (X.tmp[n+1, ] - mu[th, ]) / nvec.tmp[th]

               sigma.tmp <- sigmaSym(X = X.tmp, Y = Y.tmp, L = L,
                                               dimension = dimension, n = n + 1,
                                               nvec = nvec.tmp, B = NULL, nu=0,
                                               delta=10^(-7), prewhitened=TRUE,
                                               steps=FALSE, nmax=500)$S
             } )

      
      # delta.tmp[b] := sigma^(-1) (mu[b] - mu[th])
      delta.tmp <- t(qr.solve(sigma.tmp, t(mu.tmp) - mu.tmp[th, ]))
      # mu.th[b] := (mu.tmp[th, ] + mu.tmp[b, ]) / 2
      mu.th <- t(t(mu.tmp) + mu.tmp[th, ]) / 2

      T <- rep(0, n + 1)
      for(j in thIndex) {
        for(b in seq_len(L)[-th]) {
          T[j] <- T[j] + nvec.tmp[b] *
                          exp((X.tmp[j, ] - mu.th[b, ]) %*% delta.tmp[b, ])
        }
      }
      Tv <- T[thIndex]
      
      PV[i,th] <- sum(Tv >= tail(Tv, 1)) / nvec.tmp[th]
    }
  }
  dimnames(PV)[[2]] <- Ylevels
  return(PV)  
}

