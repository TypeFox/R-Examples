targetselection <- 
  function(data,centered=FALSE){
    if(!is.matrix(data)) data <- as.matrix(data)
    p <- nrow(data)
    N <- ncol(data)
    if(!centered){
      if(N < 4) stop("The number of columns should be greater than 3") 
    datacen <- data-rowMeans(data)
    Sigmasam <- tcrossprod(datacen)/(N-1)
    samvar <- diag(Sigmasam)    
    trSigmahat <- sum(samvar)
    Q <- sum(colSums(datacen^2)^2)/(N-1)
    trSigma2hat <- (N-1)/(N*(N-2)*(N-3))*((N-1)*(N-2)*sum(Sigmasam^2)+(trSigmahat)^2-N*Q)
    lambda1 <- (trSigmahat^2+trSigma2hat)/(N*trSigma2hat+(p-N+1)/p*trSigmahat^2)
    lambda1 <- min(lambda1,1)
    lambda2 <- (trSigmahat^2+trSigma2hat)/(N*trSigma2hat+trSigmahat^2-2*trSigmahat*(N-1)+p*(N-1))
    lambda2 <- max(0,min(lambda2,1))
    help1 <- help21 <- help22 <- help3 <- rep(0,p)
    for(i in 1:(N-1)){
      data2 <- matrix(data[,(i+1):N],p,N-i)
      help1 <- rowSums(data[,i]*data2)+help1
      help21 <- rowSums(data[,i]^3*data2)+help21
      help22 <- rowSums(data2^3*data[,i])+help22
      help3 <- rowSums(data[,i]^2*data2^2)+help3
    }
    trD21 <- 2*sum(help3)/N/(N-1)
    trD22 <- 2*(sum(help1*rowSums(data^2))-sum(help21+help22))
    trD23 <- 4*(sum(help1^2)-sum(help3)-trD22)
    trD22 <- trD22/N/(N-1)/(N-2)
    trD23 <- trD23/N/(N-1)/(N-2)/(N-3)
    trDs2 <- trD21-2*trD22+trD23   
    lambda3 <- (trSigmahat^2+trSigma2hat-2*trDs2)/(N*trSigma2hat+trSigmahat^2-(N+1)*trDs2)
    lambda3 <- max(0,min(lambda3,1))
  } else {
    if(N < 2) stop("The number of columns should be greater than 1")              
    Sigmasam <- tcrossprod(data)/N
    samvar <- diag(Sigmasam)
    trSigmahat <- sum(samvar)
    trSigma2hat <- 0
    trSigma2hat <- trDs2 <- 0
    for(i in 1:(N-1)) {
      trSigma2hat <- sum(crossprod(data[,i],data[,(i+1):N])^2) + trSigma2hat              
      trDs2 <- sum((data[,i]*data[,(i+1):N])^2) + trDs2
    }
    trSigma2hat <- 2*trSigma2hat/N/(N-1)
    trDs2 <- 2*trDs2/N/(N-1)              
    lambda1 <- (trSigmahat^2+trSigma2hat)/((N+1)*trSigma2hat+(p-N)/p*trSigmahat^2)
    lambda1 <- min(lambda1,1)
    lambda2 <- (trSigmahat^2+trSigma2hat)/((N+1)*trSigma2hat+trSigmahat^2-2*trSigmahat*N+p*N)
    lambda2 <- max(0,min(lambda2,1))
    lambda3 <- (trSigmahat^2+trSigma2hat-2*trDs2)/((N+1)*trSigma2hat+trSigmahat^2-(N+2)*trDs2)
    lambda3 <- max(0,min(lambda3,1))
  }
  cat("OPTIMAL SHRINKAGE INTENSITIES FOR THE TARGET MATRIX WITH", "\n")
  cat("Equal variances   :",round(lambda1,4),"\n")
  cat("Unit variances    :",round(lambda2,4),"\n")
  cat("Unequal variances :",round(lambda3,4),"\n")
  cat("\nSAMPLE VARIANCES", "\n")
  cat("Range   :",round(max(samvar)-min(samvar),4),"\n")
  cat("Average :",round(mean(samvar),4),"\n")
  }
