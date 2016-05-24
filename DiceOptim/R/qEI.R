qEI <- function(x,model ,plugin = NULL, type = "UK", minimization = TRUE, fastCompute = TRUE, eps = 10^(-5), envir = NULL) {
  d <- model@d
  if (!is.matrix(x)) {
    x <- matrix(x,ncol = d)
  }
  xb <- rbind(model@X,x)
  nExp <- model@n
  x <- unique(round(xb,digits = 8))
  if ((nExp+1)>length(x[,1])) {return (0)}
  x <- matrix(x[(nExp+1):length(x[,1]),],ncol=d)
  if(is.null(plugin) &&  minimization) plugin <- min(model@y)
  if(is.null(plugin) && !minimization) plugin <- max(model@y)

  if (length(x[,1]) == 1) {
    return(EI(x=x,model=model,plugin=plugin,type=type,envir=envir))
  }

  krig  <- predict(object=model, newdata=x,type=type,se.compute=FALSE, cov.compute=TRUE,checkNames = FALSE)

  sigma <- krig$cov
  mu <- krig$mean
  
  q <- length(mu)
  
  pk <- first_term <- second_term <- rep(0,times=q)
  if (!fastCompute) {
    symetric_term <- matrix(0,q,q)
    non_symetric_term <- matrix(0,q,q)
  }
  for(k in 1:q){
    #covariance matrix of the vector Z^(k)
    Sigma_k <- covZk( sigma=sigma, index=k )
    
    #mean of the vector Z^(k)
    mu_k <- mu - mu[k] 
    mu_k[k] <- -mu[k]
    if(minimization) mu_k <- -mu_k
    
    b_k <- rep(0,times=q)
    b_k[k] <- -plugin
    if(minimization) b_k <- -b_k
    pk[k] <- pmnorm(x= b_k - mu_k  ,varcov=Sigma_k,maxpts=q*200)[1]
    
    first_term[k] <- (mu[k] - plugin)*pk[k]
    if(minimization) first_term[k] <- -first_term[k]
    
    if (fastCompute) {
	    second_term[k] <- 1/eps*(pmnorm(x= b_k + eps*Sigma_k[,k]- mu_k  ,varcov=Sigma_k,maxpts=q*200)[1] - pk[k])
    } else {
        for(i in 1:q){
        	non_symetric_term[k,i] <- Sigma_k[i,k]
	       if (i>=k) {
          mik <- mu_k[i]
	         sigma_ii_k <- Sigma_k[i,i]
	            bik <- b_k[i]
	            phi_ik <- dnorm(x=bik, mean=mik , sd= sqrt(sigma_ii_k))
	  
	            #need c.i^(k) and Sigma.i^(k)
	            cik <- get_cik(b = b_k , m = mu_k , sigma = Sigma_k , i=i  )
	            sigmaik <- get_sigmaik( sigma = Sigma_k , i=i )
	            Phi_ik <- pmnorm(x = cik  ,varcov=sigmaik,maxpts=(q-1)*200)[1]
	            symetric_term[k,i] <- phi_ik*Phi_ik
	          }
          }
        }
  }
  if (!fastCompute) {
    symetric_term <- symetric_term + t(symetric_term)
    diag(symetric_term) <- 0.5 * diag(symetric_term) 
    second_term <- sum(symetric_term*non_symetric_term)
  }  
  if (!is.null(envir)) {
    assign("pk",pk,envir=envir)
    if (fastCompute == FALSE) {
      assign("symetric_term",symetric_term,envir = envir)
    }
    assign("kriging.mean", mu, envir = envir)
    assign("kriging.cov", sigma, envir = envir)
    assign("Tinv.c", krig$Tinv.c, envir = envir)
    assign("c", krig$c, envir = envir)
  }
  somFinal <- sum(first_term,second_term)
  return( somFinal )
}



get_sigmaik <- function( sigma , i ){
  
  result <- sigma
  q <- nrow(sigma)
  
  for(u in 1:q){
    for(v in 1:q){
      if((u!=i) && (v!=i)){
        result[u,v] <- sigma[u,v] - sigma[u,i]*sigma[v,i]/sigma[i,i] 
      }else{
        result[u,v] <- 0
      }
    }
  }

  result <- result[-i,-i]

  return(result)
}


get_cik <- function( b , m , sigma , i  ){
  
  sigmai <- sigma[i,] / sigma[i,i]
  
  result <- (b - m) - ( b[i] - m[i] ) * sigmai
  result <- result[-i]
  
  return(result)
  
}



covZk <- function(sigma,index){
  
  result <- sigma
  q <- nrow(sigma)
  
  result[index,index] <- sigma[index,index]
  for(i in 1:q){
    if(i!=index){
      result[index,i] <- result[i,index] <- sigma[index,index] - sigma[index,i]
    }
  }
  for(i in 1:q){
    for(j in i:q){
      if((i!=index)&&(j!=index)){
        result[i,j] <- result[j,i] <- sigma[i,j] + sigma[index,index] - sigma[index,i] - sigma[index,j]
      }
    }
  }
  
  return(result)
  
}
