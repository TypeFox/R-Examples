clustMD <-
function(X, G, CnsIndx, OrdIndx, Nnorms, MaxIter, model, store.params=FALSE, scale=FALSE, startCL="kmeans"){
    # Controls
    Y <- as.matrix(X)
    N <- nrow(Y)
    J <- ncol(Y)
    
    ### Scale continuous data
    if(scale){
      if(CnsIndx > 0)
        Y[, 1:CnsIndx] <- scale(Y[, 1:CnsIndx])
    } # if
    
    # Number of levels foon each item
    K <- apply(Y, 2, max)
    if(CnsIndx > 0) K[1:CnsIndx] <- NA
    
    # Dimension of latent space
    D <- J
    if(J > OrdIndx)
      D <- OrdIndx + sum(K[(OrdIndx+1):J] - 1)
    
    # Which dimensions correspond to each item
    if(J > OrdIndx){
      nom.ind.Z<-vector("list", J-OrdIndx)
      for(j in 1:(J-OrdIndx)){
        if(j==1){
          start <- OrdIndx + 1
        }else{
          start <- OrdIndx + sum(K[(OrdIndx+1):(OrdIndx+j-1)]-1) + 1
        }
        finish <- start + K[OrdIndx+j] - 2
        nom.ind.Z[[j]] <-c(start:finish) 
      } # j
    } # if
    
    
    ### Initial Values
    
    
    ### Estimated starting values
    ## Expected value of latent data
    Ez <- array(NA, c(N, D, G))
    for(g in 1:G) Ez[, 1:J , g] <- Y
    
    if(OrdIndx > CnsIndx){
      perc.cut <- perc.cutoffs(CnsIndx, OrdIndx, Y, N)
      zlimits <- array(NA, c(N, J, 2))
      zlimits[, 1:CnsIndx, 1] <- -Inf
      zlimits[, 1:CnsIndx, 2] <- Inf
      for(j in (CnsIndx+1):OrdIndx){
        for(k in 1:K[j]){
          zlimits[Y[,j]==k, j, 1] <- perc.cut[[j]][k]
          zlimits[Y[,j]==k, j, 2] <- perc.cut[[j]][k+1]
        }
      }
    }else{
      perc.cut <- list()
      zlimits <- array(NA, c(N, J, 2))
    }
    
    Zstart <- function(Kj, y){
      new.z <- rep( 0, (Kj - 1) )  
      if (y==1){
        new.z <- rtnorm((Kj - 1), mean=0, sd=1, upper=0)
      }else{
        new.z[-(y-1)] <- rnorm((Kj - 2), mean=0, sd=1)
        new.z[(y-1)] <- rtnorm(1, mean=0, sd=1,lower=max(new.z))
      }
      new.z
    }
    
    Zinit <- matrix(NA, N, D)
    Zinit[, 1:OrdIndx] <- Y[, 1:OrdIndx]
    if(J > OrdIndx){
      for (j in (OrdIndx+1):J){
        for(i in 1:N){
          Zinit[i, nom.ind.Z[[j-OrdIndx]]] <- Zstart(K[j], Y[i, j])
        }# i
      }# j
    }
    
    # initial clustering
    if(startCL == "kmeans"){
      if(CnsIndx > 0){
        ind <- kmeans(Y[, 1:CnsIndx], G)$cl
      }else{
        ind <- kmeans(Y, G)$cl
      }
    }else if(startCL == "hclust"){
#       temp <- hclust(dist(Y[, 1:CnsIndx]))
      temp <- hclust(dist(Y))
      ind <- cutree(temp, G)
#       print(table(ind, ind.true))
    }else if(startCL == "mclust"){
#       ind <- Mclust(Y[, 1:CnsIndx], G, model)$cl
      ind <- Mclust(Y, G, model)$cl
#       print(table(ind, ind.true))
    }else if(startCL == "random"){
      ind <- sample(1:G, N, replace=TRUE)
#       print(table(ind, ind.true))
    }
    # mixing weights
    pi.vec <- table(ind)/N
    
    # mean
    mu <- matrix(NA, D, G)
    for(g in 1:G)
      mu[, g] <- apply(Zinit[ind==g, ], 2, mean)
    
    # Covaraince
    Sigma <- array(NA, c(D, D, G))
    for(g in 1:G)
      Sigma[, , g] <- diag(D)
    
    a <- matrix(1, G, D)  
    
    ## Storage
    if(store.params==TRUE){
      ind.store <- matrix(NA, N, MaxIter)
      Ez.store <- array(NA, c(N, D, G, MaxIter))
      tau.store <- array(NA, c(N, G, MaxIter))
      mu.store <- array(NA, c(D, G, MaxIter))
      lambda.store <- array(NA, c(G, D, MaxIter))
      a.store <- array(NA, c(G, D, MaxIter))
      likeStore <- rep(NA, MaxIter)
      if(J > OrdIndx) probs.nom.store<- array(NA, c(J-OrdIndx, max(K[(OrdIndx+1):J]), G, MaxIter))
    }
    
    ### EM Loop
    for(iter in 1:MaxIter){  
      if(iter%%10==0) print(iter)
      
      # Standard normal deviates for MC approximation
      if(J > OrdIndx) norms <- mvrnorm(Nnorms, mu=rep(0, max(K[(OrdIndx+1):J])-1), Sigma=diag(max(K[(OrdIndx+1):J])-1))
      
      if(J > CnsIndx){
        # Z moments
        temp.z <- z.moments(D, G, N, CnsIndx, OrdIndx, zlimits, mu, Sigma, Y, J, K, norms, nom.ind.Z)
        Ez <- temp.z[[1]]
        S <- temp.z[[2]]
        probs.nom <- temp.z[[3]]
        S2 <- temp.z[[4]]
      }
      
      # E-step
      temp.E <- E.step(N, G, pi.vec, Y, OrdIndx, CnsIndx, D, perc.cut, mu, Sigma, Ez, J, probs.nom, K)
      tau <- temp.E[[1]]
      Elz <- temp.E[[2]]
      ind <- map(tau)
      
      # M-step
      temp.M <- M.step(tau, N, Elz, J, OrdIndx, D, G, Y, CnsIndx, S2, model, a)
      pi.vec <- temp.M[[1]]
      mu <- temp.M[[2]]
      lambda <- temp.M[[3]]
      a <- temp.M[[4]]
      Sigma <-temp.M[[5]]
      
      if(store.params==TRUE){
        ind.store[, iter] <- ind
        #       Ez.store[, , ,iter] <- Ez
        tau.store[, , iter] <- tau
        mu.store[, , iter] <- mu
        lambda.store[, , iter] <- lambda
        a.store[, , iter] <- a
        if(J > OrdIndx) probs.nom.store[, , ,iter] <- probs.nom
        likeStore[iter] <- ObsLogLikelihood(N, CnsIndx, G, Y, mu, Sigma, pi.vec, J, OrdIndx, K, perc.cut, Nnorms, zlimits, nom.ind.Z)
      }
    } # iter
    
    # approximated BIC
    obslike <- ObsLogLikelihood(N, CnsIndx, G, Y, mu, Sigma, pi.vec, J, OrdIndx, K, perc.cut, Nnorms, zlimits, nom.ind.Z)
    BIChat <- 2*obslike - npars_clustMD(model, D, G, J, OrdIndx)*log(N)
    
    if(store.params==TRUE){
      params.store.list <- list(cl.store=ind.store, tau.store=tau.store, means.store=mu.store, A.store=a.store, lambda.store=lambda.store, likelihood.store=likeStore)
      list(cl=ind, means=mu, A=a, Lambda=lambda, Sigma=Sigma, BIChat = BIChat, paramlist=params.store.list)
    }else{
      list(cl=ind, tau=tau, means=mu, A=a, Lambda=lambda, Sigma=Sigma, BIChat = BIChat)
    }
  }
