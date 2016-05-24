dhsic.test <- function(X, Y, alpha=0.05, method = "gamma", kernel = "gaussian", B = 100, pairwise = FALSE){
# Copyright (c) 2010 - 2013  Jonas Peters  [peters@stat.math.ethz.ch]
#                            Niklas Pfister [pfisteni@student.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 

  # outputs the test statistic and the critical value (according to alpha). If the test statistic is 
  # larger than the critical value, 
  # H_0 (X[[1]],...,X[[d]] are jointly independent) is rejected.

  #####################################################################################

  
  
  
  ############
  # Pairwise #
  ############
  
  if(pairwise){
    if(!missing(Y)||!is.list(X)||length(X)<=2){
      stop("pairwise only makes sense if number of variables is greater than 2")
    }
    d <- length(X)
    for(j in 1:d){
      if(is.matrix(X[[j]])==FALSE){
        X[[j]] <- as.matrix(X[[j]])
      }
    }   
    len <- nrow(X[[1]])
    
    # Case: len<2d
    if(len<2*d){
      warning("Sample size is smaller than twice the number of variables. Test is trivial.")
      test=list(statistic=0,
                crit.value = Inf,
                p.value = 1,
                time = c(GramMat=0,dHSIC=0,CritVal=0))
      return(test)
    }
    
    pValVec <- rep(0,d-1)
    statVec <- rep(0,d-1)
    critVec <- rep(0,d-1)
    ptm <- proc.time()
    for(j in (d:2)){
      resTmp <- dhsic.test(do.call("cbind", X[1:(j-1)]), X[[j]], alpha=alpha/(d-1), method = method, kernel = kernel, B = B, pairwise=FALSE)
      pValVec[d-j+1] <- resTmp$p.value
      statVec[d-j+1] <- resTmp$statistic
      critVec[d-j+1] <- resTmp$crit.value
    }
    ind <- which.min(pValVec)
    stat <- statVec[ind]
    critical_value <- critVec[ind]
    # bonferroni correction for p-value
    p_value = min(pValVec[ind]*(d-1),1)
    timeCritVal <- as.numeric((proc.time()-ptm)[1])
    timeGramMat <- as.numeric(resTmp$time[1])
    timeHSIC <- as.numeric(resTmp$time[2])
  }

  ###################
  # d-variable HSIC #
  ###################

  else{
    
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
      warning("Sample size is smaller than twice the number of variables. Test is trivial.")
      test=list(statistic=0,
                crit.value = Inf,
                p.value = 1,
                time = c(GramMat=0,dHSIC=0,CritVal=0))
      return(test)
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
      cols <- ncol(x)
      if(cols>1){
        temp <- matrix(rep(x[,1],len),ncol=len)
        KX <- ((temp-t(temp))==0)*1 
        for(i in 2:cols){
          temp <- matrix(rep(x[,i],len),ncol=len)
          KX <- KX*(((temp-t(temp))==0)*1)
        }
      }
      else{
        temp <- matrix(rep(x,len),ncol=len)
        KX <- ((temp-t(temp))==0)*1 
      }
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
    bandwidth <- vector("numeric",d)
    ptm <- proc.time()
    for(j in 1:d){
      if(kernel[j]=="gaussian"){
        bandwidth[j] <- median_bandwidth(X[[j]])
        K[[j]] <- gaussian_grammat(X[[j]],bandwidth[j])
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
    # Main Computation
    ###
    
    ### Compute dHSIC
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
    
    ### Using permutation test for dHSIC
    ptm <- proc.time()
    if(method == "permutation"){
      # Define permutation function
      dhsic_perm_fun <- function(i){
        term1 <- K[[1]]
        term2 <- 1/len^2*sum(K[[1]])
        term3 <- 2/len^2*colSums(K[[1]])
        for (j in 2:d){
          X_perm <- matrix(X[[j]][sample(len,replace=FALSE),],nrow=len)
          if(kernel[j]=="gaussian"){
            K_perm <- gaussian_grammat(X_perm,bandwidth[j])
          }
          else if(kernel[j]=="discrete"){
            K_perm <- discrete_grammat(X_perm)
          }
          else{
            K_perm <- custom_grammat(X_perm,kernel[j])
          }
          term1 <- term1*K_perm
          term2 <- 1/len^2*term2*sum(K_perm)
          term3 <- 1/len*term3*colSums(K_perm)
        }
        term1 <- sum(term1)
        term3 <- sum(term3)
        return(1/len^2*term1+term2-term3)
      }
      # Perform performing
      dHSIC_perm <- sapply(1:B,dhsic_perm_fun)
      
      # Compute critical value and p-value
      sortdHSIC <- sort(len*dHSIC_perm)
      Bind <- sum(len*dHSIC==sortdHSIC)+ceiling((1-alpha)*(B+1))
      if(Bind<=B){
        critical_value <- sortdHSIC[Bind]
      }
      else{
        critical_value <- Inf
      }
      p_value <- (sum(dHSIC_perm>=dHSIC)+1)/(B+1)
    }
  
    ### Using bootstrap test for dHSIC
    if (method == "bootstrap"){
      # Define bootstrap function
      dhsic_boot_fun <- function(i){
        term1 <- K[[1]]
        term2 <- 1/len^2*sum(K[[1]])
        term3 <- 2/len^2*colSums(K[[1]])
        for (j in 2:d){
          X_boot <- matrix(X[[j]][sample(len,replace=TRUE),],nrow=len)
          if(kernel[j]=="gaussian"){
            K_boot <- gaussian_grammat(X_boot,bandwidth[j])
          }
          else if(kernel[j]=="discrete"){
            K_boot <- discrete_grammat(X_boot)
          }
          else{
            K_boot <- custom_grammat(X_boot,kernel[j])
          }
          term1 <- term1*K_boot
          term2 <- 1/len^2*term2*sum(K_boot)
          term3 <- 1/len*term3*colSums(K_boot)
        }
        term1 <- sum(term1)
        term3 <- sum(term3)
        return(1/len^2*term1+term2-term3)
      }
      # Perform bootstrapping
      dHSIC_boot <- sapply(1:B,dhsic_boot_fun)
      
      # Compute critical value and p-value
      sortdHSIC <- sort(len*dHSIC_boot)
      Bind <- sum(len*dHSIC==sortdHSIC)+ceiling((1-alpha)*(B+1))
      if(Bind<=B){
        critical_value <- sortdHSIC[Bind]
      }
      else{
        critical_value <- Inf
      }
      p_value <- (sum(dHSIC_boot>=dHSIC)+1)/(B+1)
    }
    
    ### Gamma approximation based test for dHSIC
    else if (method == "gamma"){
      # estimators
      est.a <- vector("numeric", d)
      est.b <- vector("numeric", d)
      est.c <- vector("numeric", d)
      for (j in 1:d){
        est.a[j] <- 1/(len^2)*sum(K[[j]])
        est.b[j] <- 1/(len^2)*sum(K[[j]]^2)
        est.c[j] <- 1/len^3*sum(rowSums(K[[j]])^2)
      }
      prod.a <- prod(est.a)
      prod.b <- prod(est.b)
      prod.c <- prod(est.c)
      oneoutprod.a <- vector("numeric",d)
      oneoutprod.b <- vector("numeric",d)
      oneoutprod.c <- vector("numeric",d)
      for(j in 1:d){
        oneoutprod.a[j] <- prod.a/est.a[j]
        oneoutprod.b[j] <- prod.b/est.b[j]
        oneoutprod.c[j] <- prod.c/est.c[j]
      }
      est.d <- est.a^2
      prod.d <- prod.a^2
      oneoutprod.d <- oneoutprod.a^2
      # exp-estimator
      exp.est <- (1-sum(oneoutprod.a)+(d-1)*prod.a)/len
      # var-estimtor
      term1 <- prod.b
      term2 <- (d-1)^2*prod.d
      term3 <- 2*(d-1)*prod.c
      term4 <- 0
      term5 <- 0
      term6 <- 0
      term7 <- 0
      for (r in 1:(d-1)){
        term4 <- term4+est.b[r]*oneoutprod.d[r]
        term5 <- term5+est.b[r]*oneoutprod.c[r]
        term6 <- term6+est.c[r]*oneoutprod.d[r]
        for (s in (r+1):d){
          term7 <- term7+2*est.c[r]*est.c[s]*oneoutprod.d[r]/est.d[s]
        }
      }
      term4 <- term4+est.b[d]*oneoutprod.d[d]
      term5 <- -2*(term5+est.b[d]*oneoutprod.c[d])
      term6 <- -2*(d-1)*(term6+est.c[d]*oneoutprod.d[d])
      
      factor1 <- len-2*d
      factor2 <- len*(len-1)*(len-2)
      for (j in 1:(2*d-3)){
        factor1 <- factor1*(len-2*d-j)
        factor2 <- factor2*(len-2-j)
      }
      var.est=2*factor1/factor2*(term1+term2+term3+term4+term5+term6+term7)
      # Compute critical value and p-value
      a <- (exp.est^2)/var.est
      b <- len*var.est/exp.est
      critical_value <- qgamma(1-alpha,shape=a,scale=b)
      p_value <- pgamma(len*dHSIC,shape=a,scale=b, lower.tail=FALSE)
    }
    
    
    ### Eigenvalue based test for dHSIC
    else if (method == "eigenvalue"){
      ## Use eigenvalue approximation to calculate the critical value
      # calculate H2
      est1 <- vector("numeric",d)
      est2 <- matrix(NA,len,d)
      for (j in 1:d){
        est1[j] <- 1/(len*(len-1))*sum(K[[j]]-diag(diag(K[[j]])))
        est2[,j] <- 1/len*rowSums(K[[j]])
      }
      est1_prod <- prod(est1)
      est2_prod <- apply(est2,1,prod)
      
      a1 <- matrix(1,len,len)
      a2 <- (d-1)^2*est1_prod
      a3 <- (d-1)*est2_prod
      a5 <- matrix(0,len,len)
      a6 <- matrix(0,len,len)
      a8 <- matrix(0,len,len)
      a9 <- matrix(0,len,1)
      for (j in 1:d){
        a1 <- a1*K[[j]]
        a5 <- a5+K[[j]]*est1_prod/est1[j]
        a6 <- a6+K[[j]]*matrix(rep(est2_prod/est2[,j],len),len,len)
        a9 <- a9+est2[,j,drop=FALSE]*est1_prod/est1[j]
        i <- j+1
        while ( i<=d ){
          a8 <- a8+est2[,j,drop=FALSE]%*%t(est2[,i,drop=FALSE])*est1_prod/(est1[j]*est1[i])+est2[,i,drop=FALSE]%*%t(est2[,j,drop=FALSE])*est1_prod/(est1[j]*est1[i])
          i <- i+1
        }
      }
      a3 <- matrix(rep(a3,len),len,len)
      a4 <- t(a3)
      a6 <- -a6
      a7 <- t(a6)
      a8 <- a8
      a9 <- matrix(rep((1-d)*a9,len),len,len)
      a10 <- t(a9)
      
      H2 <- 1/(d*(2*d-1))*(a1+a2+a3+a4+a5+a6+a7+a8+a9+a10)
      
      eigenvalues <- eigen(H2,only.values = TRUE,symmetric = TRUE)$values/len
      
      # simulate the sum of chi-squared distribution
      M <- 5000
      Z <- rnorm(M*len)^2
      chi_dist <- d*(2*d-1)*colSums(matrix(Z*eigenvalues,len,M))
      critical_value <- as.numeric(quantile(chi_dist,1-alpha))
      p_value <- 1-ecdf(chi_dist)(dHSIC*len)
    }
    timeCritVal <- as.numeric((proc.time()-ptm)[1])
    stat <- dHSIC*len
  }
  
  ###
  # Collect result
  ###
  test=list(statistic=stat,
            crit.value = critical_value,
            p.value = p_value,
            time = c(GramMat=timeGramMat,dHSIC=timeHSIC,CritVal=timeCritVal))
  
  return(test)
}
