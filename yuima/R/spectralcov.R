
## Local Method of Moments estimator

lmm <- function(x, block = 20, freq = 50, freq.p = 10, K = 4, interval = c(0, 1),
                Sigma.p = NULL, noise.var = "AMZ", samp.adj = "direct", psd = TRUE){
  
  Data <- get.zoo.data(x)
  
  d.size <- length(Data)
  N <- sapply(Data,"length") - 1
  n <- min(N)
  
  interval <- as.numeric(interval)
  freq.p <- min(freq.p,freq)
  
  if(d.size > 1){
    
    # Constructing the stpectral statistics and noise level estimates
    Spec <- array(0,dim=c(d.size,block,freq))
    H <- matrix(0,d.size,block)
    cH <- array(0,dim=c(d.size^2,block,freq))
    zmat <- diag(d.size^2)
    
    phi.br <- (block*pi*(1:freq))^2
    
    for(p in 1:d.size){
      
      dY <- diff(as.numeric(Data[[p]]))
      Time <- (as.numeric(time(Data[[p]])) - interval[1])/(interval[2] - interval[1])
      
      if((min(Time) < 0) || (max(Time) > 1)){
        print("lmm:invalid interval")
        return(NULL)
      }
      
      Time2 <- (Time[1:N[p]]+Time[-1])/2
      
      k.idx <- findInterval(Time2,vec=seq(0,1,by=1/block))
      
      # Noise variance
      if(is.numeric(noise.var)){
        sq.eta <- noise.var[p]
      }else{
        sq.eta <- switch(noise.var,
                         "BR" = mean(dY^2)/2, # Bandi and Russel (2006)
                         "O" = -sum(dY[-N[p]] * dY[-1])/N[p], # Oomen (2006)
                         "AMZ" = {tmp <- arima0(dY,order=c(0,0,1),include.mean=FALSE);
                         -tmp$coef * tmp$sigma2} # Ait-Sahalia et al. (2005)
                         )
      }
      
      # Effect of irregularity
      if(is.numeric(samp.adj)){
        EOI <- samp.adj[p, ]
      }else{
        EOI <- switch(samp.adj,
                      "QVT" = block * rowsum(diff(Time)^2,group = findInterval(Time[-1],vec=seq(0,1,by=1/block),
                                                                               rightmost.closed = TRUE)),
                      "direct" = block * rowsum((diff(Time,lag=2)/2)^2, group = k.idx[-N[p]]))
      }
      
      Spec[p,,] <- sqrt(2*block)*
        rowsum(sin(tcrossprod(pi*(block*Time2-(k.idx-1)),1:freq))*dY,
               group = k.idx)
      H[p, ] <- sq.eta * EOI
      cH[d.size*(p-1)+p,,] <- tcrossprod(H[p,],phi.br)
      
      for(q in 1:d.size){
        
        tmp.mat <- diag(0, d.size)
        tmp.mat[p, q] <- 1
        zmat <- zmat + tmp.mat %x% t(tmp.mat)
        
      }
      
    }
    
    #summand <- apply(Spec,c(2,3),FUN="tcrossprod") - cH
    summand <- array(.C("krprod",
                        as.double(Spec),
                        as.integer(d.size),
                        as.integer(block*freq),
                        result=double(d.size^2*block*freq))$result,
                     dim=c(d.size^2,block,freq)) - cH
    
    if(is.null(Sigma.p)){
      # Pilot estimation of Sigma
      # Using K adjacent blocks following Bibinger et al. (2014a)
      #Sigma.p <- rollmean(t(rowMeans(array(summand[,,1:freq.p],
      #                                      dim=c(d.size^2,block,freq.p)),
      #                                dims = 2)),k = K, fill = "e")
      # Recent version of Bibinger et al. (2014b)
      Sigma.p <- rollapply(t(rowMeans(array(summand[,,1:freq.p],
                                            dim=c(d.size^2,block,freq.p)),
                                      dims = 2)), width = 2*K + 1, 
                           FUN="mean", partial = TRUE)
    }
    
    # Local method of moments
    Sigma <- matrix(0,d.size^2,block)
    inv.Ik <- array(0,dim=c(d.size^2,d.size^2,block))
    
    for(k in 1:block){
      
      tmp <- double(d.size^2)
      H.inv <- diag(1/H[,k])
      tmp.eigen <- eigen(matrix(Sigma.p[k,], d.size, d.size) %*% H.inv,
                         symmetric=FALSE)
      inv.V <- solve(tmp.eigen$vectors)
      iHV <- H.inv %*% tmp.eigen$vectors
      iHVxiHV <- iHV %x% iHV
      inv.VxV <- inv.V %x% inv.V
      #Lambda <- apply(1/outer(tmp.eigen$values, phi.br, FUN="+"),
      #                2,FUN="tcrossprod")
      Lambda <- matrix(.C("krprod",
                          as.double(1/outer(tmp.eigen$values, phi.br, FUN="+")),
                          as.integer(d.size),
                          as.integer(freq),
                          result=double(d.size^2*freq))$result,d.size^2,freq)
      #Lambda <- KhatriRao(1/outer(tmp.eigen$values, phi.br, FUN="+"))
      inv.Ik[,,k] <- solve(iHVxiHV %*% diag(rowSums(Lambda)) %*% inv.VxV)
      Sigma[,k] <- inv.Ik[,,k] %*% iHVxiHV %*% 
        rowSums(Lambda * inv.VxV %*% summand[,k,])
      
    }
    
    cmat <- matrix(rowMeans(Sigma), d.size, d.size)
    vcov <- rowMeans(inv.Ik, dims = 2) %*% zmat/block
    
    if(psd){
      r <- eigen(cmat, symmetric = TRUE)
      cmat <- r$vectors %*% abs(diag(r$values)) %*% t(r$vectors)
      r <- eigen(vcov, symmetric = TRUE)
      vcov <- r$vectors %*% abs(diag(r$values)) %*% t(r$vectors)
    }
    
    result <- list(covmat = cmat, vcov = vcov, Sigma.p = Sigma.p)
    class(result) <- "yuima.specv"
    
  }else{
    
    # Constructing the stpectral statistics and noise level estimates
    dY <- diff(as.numeric(Data[[1]]))
    Time <- (as.numeric(time(Data[[1]])) - interval[1])/(interval[2] - interval[1])
    
    if((min(Time) < 0) || (max(Time) > 1)){
      print("lmm:invalid interval")
      return(NULL)
    }
    
    Time2 <- (Time[1:N]+Time[-1])/2
    
    k.idx <- findInterval(Time2,vec=seq(0,1,by=1/block))
    
    # Noise variance
    if(is.numeric(noise.var)){
      sq.eta <- noise.var
    }else{
      sq.eta <- switch(noise.var,
                       "BR" = mean(dY^2)/2, # Bandi and Russel (2006)
                       "O" = -sum(dY[-N] * dY[-1])/N, # Oomen (2006)
                       "AMZ" = {tmp <- arima0(dY,order=c(0,0,1),include.mean=FALSE);
                       -tmp$coef * tmp$sigma2} # Ait-Sahalia et al. (2005)
      )
    }
    
    phi.br <- (block*pi*(1:freq))^2
    
    # Effect of irregularity
    if(is.numeric(samp.adj)){
      EOI <- samp.adj
    }else{
      EOI <- switch(samp.adj,
                    "QVT" = block * rowsum(diff(Time)^2,group = findInterval(Time[-1],vec=seq(0,1,by=1/block),
                                                                             rightmost.closed = TRUE)),
                    "direct" = block * rowsum((diff(Time,lag=2)/2)^2, group = k.idx[-N]))
    }
    
    Spec <- sqrt(2*block)*
      rowsum(sin(tcrossprod(pi*(block*Time2-(k.idx-1)),1:freq))*dY,
             group = k.idx)
    H <- sq.eta * EOI
    cH <- tcrossprod(H,phi.br)
    
    summand <- Spec^2 - cH
    
    if(is.null(Sigma.p)){
      # Pilot estimation of Sigma
      # Using K adjacent blocks following Bibinger et al. (2014a)
      #Sigma.p <- rollmean(rowMeans(matrix(summand[,1:freq.p],block,freq.p)),
      #                    k = K, fill = "e")
      # Recent version of Bibinger et al. (2014b)
      Sigma.p <- rollapply(rowMeans(matrix(summand[,1:freq.p], block,freq.p)), 
                           width = 2*K + 1, FUN="mean", partial = TRUE)
    }
    
    Ijk <- 1/(Sigma.p + cH)^2
    inv.Ik <- 1/rowSums(Ijk)
    
    if(psd){
      cmat <- sum(Ijk*inv.Ik*summand)/block
    }else{
      cmat <- abs(sum(Ijk*inv.Ik*summand)/block)
    }
    
    result <- list(covmat = as.matrix(cmat), vcov = as.matrix(2 * mean(inv.Ik)/block), Sigma.p = Sigma.p)
    class(result) <- "yuima.specv"
    
  }
  
  return(result)
}


# print method for yuima.specv-class
print.yuima.specv <- function(x, ...){
  
  cat("Estimated covariance matrix\n")
  print(x$covmat, ...)
  cat("Standard Error\n")
  print(matrix(sqrt(diag(x$vcov)), ncol = ncol(x$covmat)), ...)
  
}

