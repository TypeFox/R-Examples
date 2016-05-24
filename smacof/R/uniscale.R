## Performs unidimensional scaling based on all possible dissimilarity permutations
## It returns the configuration based on the stress minimum. 


uniscale <- function(delta, weightmat = NULL) {
  delta <- as.matrix(delta)
  if (is.null(weightmat)) {
    w <- 1-diag(dim(delta)[1])
  } else {
    w <- as.matrix(weightmat)
  }
  n <- dim(delta)[1]
  nn <- n*(n-1)/2
  
  diss <- delta
  delta <- as.matrix(normDissN(as.dist(delta), as.dist(w), 1))
  
  m <- 0
  k <- 0
  fmin <- Inf
  x <- 1:n
  v <- as.matrix(solve((diag(rowSums(w))-w)+(1/n))-(1/n))
  
  
  repeat{
    k <- k+1
 	  s <- sign(outer(x,x,"-"))
	  t <- as.vector(v%*%rowSums(delta*w*s))
	  if (are.monotone(x,t)) {
	 	 m <- m+1
		 d <- abs(outer(t,t,"-"))   
		 f <- sum(as.dist(w*(delta-d)^2))/nn
     if (f < fmin) {
		 	 fmin <- f
		 	 xmin <- t
		 }
		}
	if (prod(x==(n:1))==1) break ##return(list(xmin = xmin, fmin = fmin, m = m, k = k))
	x <- next.perm(x)
	}
  
  confdiss <- dist(xmin)          ## configuration distances
  stress <- sqrt(fmin)         ## stress-1 normalization (/2 because of full matrix instead of dist in loop)
  names(xmin) <- colnames(diss)
  
  res <- list(delta = as.dist(diss), conf = xmin, confdiss = as.dist(confdiss), stress = stress, weightmat = as.dist(w), npermtot = k, 
              npermscale = m, nobj = n, call = match.call())
  class(res) <- "uniscale"
  return(res)
}
	