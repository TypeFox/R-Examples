InitParameter <- 
  function(data, 
           nbClust, 
           n.start = 250, 
           small.pen = 0.5)
  { 
    data <- as.matrix(scale(data, TRUE, FALSE))
    n <- as.integer(dim(data)[1])
    p <- as.integer(dim(data)[2])
    n.start <- as.integer(n.start)
    nbClust <- as.integer(nbClust)
    small.pen <- as.double(small.pen)
    
    
    Mu <- matrix(0, p, nbClust)
    clust <- kmeans(data, nbClust, iter.max = 1e3, nstart = n.start, algorithm = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"))
    memb <- clust$cluster
    
    S <- array(0, dim = c(p, p, nbClust))
    for(k in 1:nbClust)
    {
      Mu[,k] <- colMeans(data[memb == k, , drop = FALSE])
      S[,,k]  <- cov(data[memb == k, , drop = FALSE])
    }
    W <- Wi <- array(0, dim = c(p, p, nbClust))
    for(k in 1:nbClust)
    {
      gg <- glasso(S[,,k], rho = small.pen,  thr=1e-4, maxit=1e4,  approx=FALSE, penalize.diagonal=FALSE, start=c("cold","warm"), w.init = NULL, wi.init=NULL, trace=FALSE) 
      W[,,k] <- gg$w
      Wi[,,k] <- gg$wi
    }
    
    prop <- clust$size/n
    P <- list()
    P$X <-  data
    P$prop <-  prop
    P$Mu <-  Mu
    P$CovarianceMatrix <- W
    P$PrecisionMatrix <-  Wi
    P$ProbCond <- matrix(0, nrow(data), nbClust)
    return(P)
  }
