kGmedian <- function (X, ncenters=2, gamma=1, alpha=0.75, nstart = 10, nstartkmeans=10) 
{
  ### 
  X <- as.matrix(X)
  p <- ncol(X)
  n <- nrow(X)
  
  ### Initialization  
  if (is.matrix(ncenters)==FALSE) {
    k <- ncenters
    centers = kmeans(X, ncenters,nstart=nstartkmeans,algorithm="MacQueen")$centers
    
  } 
  else {
    k <- nrow(ncenters)
    centers <- ncenters
  }
  Z <- stoKmed_rcpp(X, X, centers, gamma=gamma, alpha = alpha)
  best <- sum(Z$wss)
  
  if (nstart >= 2) {
    for (i in 2:nstart) {
      #      ind.cent = sample(c(k:(m+k)),k)
      #			centers = x[ind.cent,]		
      ind.init = sample(c(1:n),k)
      centers = X[ind.init, ]
      x0 = X[-ind.init,]
      ZZ <-  stoKmed_rcpp(x0, X, centers, gamma=gamma, alpha = alpha)
      if ((z <- sum(ZZ$wss)) < best) {
        Z = ZZ
        best = z
      }
    }
  }
  centers = matrix(Z$centers, k)
  dimnames(centers) = list(1L:k, dimnames(X)[[2L]])
  cluster = Z$cl
  #if (!is.null(rn <- rownames(x))) 
  #        names(cluster) <- rn
  out <- list(cluster = cluster, centers = centers, withinsrs = Z$wss, 
              size = Z$nc)
  # class(out) <- "kmeans"*/
  return(out)
}
