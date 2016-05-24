dhsic <- function(X, Y, kernel = "gaussian"){
# Copyright (c) 2010 - 2013  Jonas Peters  [peters@stat.math.ethz.ch]
#                            Niklas Pfister [pfisteni@student.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 

  # outputs the empirical pHSIC

  ###
  # Prerequisites
  ###
  
  # Check if Y has been passed as argument (implies two variable case)
  if(!missing(Y)){
     X <- list(X,Y)
  }
      
  # Set d, convert to matricies and set len
  d <- length(X)
  for(j in 1:d){
    if(is.matrix(X[[j]])==FALSE){
      X[[j]] <- as.matrix(X[[j]])
    }
  }     
  len <- nrow(X[[1]])
  
  # Case: len<2d
  if(len<2*d){
    warning("Sample size is smaller than twice the number of variables. dHSIC is trivial.")
    result=list(dHSIC = 0,
                time = c(GramMat=0,HSIC=0))
    return(result)
  }
  
  
  # Check if enough kernels where set given the number of variables (else only used first kernel)
  if(length(kernel)<d){
    kernel <- rep(kernel[1],d)
  }
  
  # Define median heuristic bandwidth function
  median_bandwidth <- function(x){
    xnorm <- as.matrix(dist(x,method="euclidean",diag=TRUE,upper=TRUE))
    if(len>1000){
      sam <- sample(1:len,1000)
      xhilf <- xnorm[sam,sam]
    }
    else{
      xhilf <- xnorm
    }
    bandwidth <- sqrt(0.5*median(xhilf[lower.tri(xhilf,diag=FALSE)]^2))
    if(bandwidth==0){
      bandwidth <- 0.001
    }
    return(bandwidth)
  }
  # Define Gaussian Gram-matrix
  gaussian_grammat <- function(x,bandwidth){
    xnorm <- as.matrix(dist(x,method="euclidean",diag=TRUE,upper=TRUE))
    xnorm <- xnorm^2
    KX <- exp(-xnorm/(2*bandwidth^2))
    return(KX)
  }
  # Define discrete kernel Gram-matrix
  discrete_grammat <- function(x){
    temp <- matrix(rep(x,len),ncol=len)
    KX <- ((temp-t(temp))==0)*1
    return(KX)
  }
  
  # Define custom kernel Gram-matrix
  custom_grammat <- function(x,fun){
    KX <- matrix(nrow=len,ncol=len)
    for (i in 1:len){
      for (j in i:len){
        KX[i,j] <- match.fun(fun)(x[i,],x[j,])
        KX[j,i] <- KX[i,j]
      }
    }
    return(KX)
  }
  

  ###
  # Compute Gram-matrix (using specified kernels)
  ###
  K <- vector("list", d)
  ptm <- proc.time()
  for(j in 1:d){
    if(kernel[j]=="gaussian"){
      bandwidth <- median_bandwidth(X[[j]])
      K[[j]] <- gaussian_grammat(X[[j]],bandwidth)
    }
    else if(kernel[j]=="discrete"){
      K[[j]] <- discrete_grammat(X[[j]])
    }
    else{
      K[[j]] <- custom_grammat(X[[j]],kernel[j])
    }
  }
  timeGramMat <- as.numeric((proc.time() - ptm)[1])

  ###
  # Compute dHSIC
  ###
  ptm <- proc.time()
  term1 <- 1
  term2 <- 1
  term3 <- 2/len
  for (j in 1:d){
    term1 <- term1*K[[j]]
    term2 <- 1/len^2*term2*sum(K[[j]])
    term3 <- 1/len*term3*colSums(K[[j]])
  }
  term1 <- sum(term1)
  term3 <- sum(term3)
  dHSIC=1/len^2*term1+term2-term3
  timeHSIC <- as.numeric((proc.time() - ptm)[1])

  ###
  # Collect result
  ###
  result=list(dHSIC = dHSIC,
              time = c(GramMat=timeGramMat,HSIC=timeHSIC))
  
  return(result)

}
