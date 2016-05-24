remlvc <-
  function(yty,Xty,Zty,XtX,ZtZ,XtZ,ndf,rdm,tau=rep(1,length(rdm)),imx=100,tol=10^-5,alg=c("FS","NR","EM")){
    ###### REML Estimation of Variance Components (via Fisher Scoring, Newton-Raphson, or EM)
    ###### Nathaniel E. Helwig (helwig@umn.edu)
    ###### Last modified: August 9, 2015
    
    ### Inputs: 
    # yty: crossprod(y)
    # Xty: crossprod(X,y)
    # Zty: crossprod(Z,y)
    # XtX: crossprod(X)
    # ZtZ: diag(crossprod(Z)) or crossprod(Z)
    # XtZ: crossprod(X,Z)
    # ndf: nrow(X) - ncol(X)
    # rdm: r by 1 random effects index dimension vector
    # tau: initial estimate of variance components
    # imx: maximum number of iterations
    # tol: convergence tolerance
    # alg: estimation algorithm
    
    ### Notes: 
    #   y = X%*%alpha + Z%*%beta + e   where
    #
    #   y is n by 1 response vector
    #   X is n by p fixed effect design matrix
    #     note: df is typically ncol(X)
    #   alpha is p by 1 fixed effect coefficients
    #   Z is n by q random effect design matrix 
    #     note: Z = cbind(Z.1, Z.2, ... , Z.r) 
    #     where Z.k is n by q.k
    #   beta is q by 1 random effect coefficients
    #     note: beta = c(beta.1, beta.2, ... ,beta.r)
    #     where beta.k ~ N(0,sig*tau.k*diag(q.k))
    #   rdm = c(q.1, q.2, ... , q.r) and q = sum_{k=1}^{r} q.k
    #   tau = c(tau.1, tau.2, ... , tau.r)
    #   e is n by 1 residual vector
    #     note: e ~ N(0,sig*diag(n))
    
    #browser(1>0)
    
    ### get initial info
    nre <- length(rdm)
    if(length(tau)!=nre){stop("Inputs 'rdm' and 'tau' must have same length.")}
    cdm <- c(0,cumsum(rdm))
    XtXi <- pinvsm(XtX)
    XtXiZ <- XtXi%*%XtZ
    ZtZX <- (-1)*crossprod(XtZ,XtXiZ)
    if(is.matrix(ZtZ)) {ZtZX <- ZtZX+ZtZ} else {diag(ZtZX) <- diag(ZtZX)+ZtZ}
    ZtyX <- Zty-crossprod(XtZ,XtXi)%*%Xty
    
    ### initialize Dmat, sig, and n2LL
    xeig <- eigen(ZtZX+diag(rep(1/tau,times=rdm)),symmetric=TRUE)
    nze <- sum(xeig$val>xeig$val[1]*.Machine$double.eps)
    Dmat <- xeig$vec[,1:nze]%*%((xeig$val[1:nze]^-1)*t(xeig$vec[,1:nze]))
    Bmat <- XtXiZ%*%Dmat
    alpha <- (XtXi+Bmat%*%t(XtXiZ))%*%Xty-Bmat%*%Zty
    beta <- Dmat%*%ZtyX
    sig <- (yty-t(c(Xty,Zty))%*%c(alpha,beta))/ndf
    if(sig<=0) sig <- 10^-3
    n2LLold <- sum(log(xeig$val[1:nze])) + sum(log(tau)*rdm) + ndf*log(sig)
    
    ### miscellanous initializations
    alg <- alg[1]
    vtol <- 1
    iter <- 0
    ival <- TRUE
    
    ### which algorithm
    if(alg=="FS"){
      # Fisher Scoring
      
      # initialize score and information
      gvec <- rep(0,nre)
      hmat <- matrix(0,nre,nre)
      for(j in 1:nre){
        idx <- (1+cdm[j]):cdm[j+1]
        trval <- sum(diag(Dmat[idx,idx]))
        gval <- (rdm[j]/tau[j])-(trval/(tau[j]^2))
        gvec[j] <- (-1)*gval+crossprod(beta[idx])/((tau[j]^2)*sig)
        hmat[j,j] <- (rdm[j]/(tau[j]^2))-2*(trval/(tau[j]^3))+sum(Dmat[idx,idx]^2)/(tau[j]^4)
        if(j>1){
          for(k in 1:(j-1)){
            kidx <- (1+cdm[k]):cdm[k+1]
            hmat[j,k] <- hmat[k,j] <- sum(Dmat[idx,kidx]^2)/((tau[j]^2)*(tau[k]^2))
          }
        }
      }
      
      # iterative update
      while(ival) {
        
        # update tau parameter estimates
        tau <- tau + pinvsm(hmat)%*%gvec
        if(any(tau<=0)) tau[tau<=0] <- 10^-3
        
        # update Dmat, sig, and n2LL
        xeig <- eigen(ZtZX+diag(rep(1/tau,times=rdm)),symmetric=TRUE)
        nze <- sum(xeig$val>xeig$val[1]*.Machine$double.eps)
        Dmat <- xeig$vec[,1:nze]%*%((xeig$val[1:nze]^-1)*t(xeig$vec[,1:nze]))
        Bmat <- XtXiZ%*%Dmat
        alpha <- (XtXi+Bmat%*%t(XtXiZ))%*%Xty-Bmat%*%Zty
        beta <- Dmat%*%ZtyX
        sig <- (yty-t(c(Xty,Zty))%*%c(alpha,beta))/ndf
        if(sig<=0) sig <- 10^-3
        n2LLnew <- sum(log(xeig$val[1:nze])) + sum(log(tau)*rdm) + ndf*log(sig)
        
        # check for convergence (and update score and information)
        vtol <- abs((n2LLold-n2LLnew)/n2LLold)
        iter <- iter + 1L
        if(vtol>tol && iter<imx){
          n2LLold <- n2LLnew
          for(j in 1:nre){
            idx <- (1+cdm[j]):cdm[j+1]
            trval <- sum(diag(Dmat[idx,idx]))
            gval <- (rdm[j]/tau[j])-(trval/(tau[j]^2))
            gvec[j] <- (-1)*gval+crossprod(beta[idx])/((tau[j]^2)*sig)
            hmat[j,j] <- (rdm[j]/(tau[j]^2))-2*(trval/(tau[j]^3))+sum(Dmat[idx,idx]^2)/(tau[j]^4)
            if(j>1){
              for(k in 1:(j-1)){
                kidx <- (1+cdm[k]):cdm[k+1]
                hmat[j,k] <- hmat[k,j] <- sum(Dmat[idx,kidx]^2)/((tau[j]^2)*(tau[k]^2))
              }
            }
          }
        } else {
          ival <- FALSE
        }
        
      } # end while(ival)
      
    } else if(alg=="NR"){
      # Newton-Raphson
      
      # initialize score and information
      gvec <- rep(0,nre)
      hmat <- matrix(0,nre,nre)
      for(j in 1:nre){
        idx <- (1+cdm[j]):cdm[j+1]
        trval <- sum(diag(Dmat[idx,idx]))
        gval <- (rdm[j]/tau[j])-(trval/(tau[j]^2))
        gvec[j] <- (-1)*gval+crossprod(beta[idx])/((tau[j]^2)*sig)
        hmat[j,j] <- 2*(trval/(tau[j]^3)) - (rdm[j]/(tau[j]^2)) - sum(Dmat[idx,idx]^2)/(tau[j]^4)
        hmat[j,j] <- hmat[j,j] + 2*crossprod(beta[idx])/((tau[j]^3)*sig) - 2*crossprod(beta[idx],Dmat[idx,idx]%*%beta[idx])/((tau[j]^4)*sig)
        if(j>1){
          for(k in 1:(j-1)){
            kidx <- (1+cdm[k]):cdm[k+1]
            hmat[j,k] <- hmat[k,j] <- (-1)*( sum(Dmat[idx,kidx]^2) + 2*crossprod(beta[idx],Dmat[idx,kidx]%*%beta[kidx])/sig ) / ((tau[j]^2)*(tau[k]^2))
          }
        }
      }
      
      # iterative update
      while(ival) {
        
        # update tau parameter estimates
        tau <- tau + pinvsm(hmat)%*%gvec
        if(any(tau<=0)) tau[tau<=0] <- 10^-3
        
        # update Dmat, sig, and n2LL
        xeig <- eigen(ZtZX+diag(rep(1/tau,times=rdm)),symmetric=TRUE)
        nze <- sum(xeig$val>xeig$val[1]*.Machine$double.eps)
        Dmat <- xeig$vec[,1:nze]%*%((xeig$val[1:nze]^-1)*t(xeig$vec[,1:nze]))
        Bmat <- XtXiZ%*%Dmat
        alpha <- (XtXi+Bmat%*%t(XtXiZ))%*%Xty-Bmat%*%Zty
        beta <- Dmat%*%ZtyX
        sig <- (yty-t(c(Xty,Zty))%*%c(alpha,beta))/ndf
        if(sig<=0) sig <- 10^-3
        n2LLnew <- sum(log(xeig$val[1:nze])) + sum(log(tau)*rdm) + ndf*log(sig)
        
        # check for convergence (and update score and information)
        vtol <- abs((n2LLold-n2LLnew)/n2LLold)
        iter <- iter + 1L
        if(vtol>tol && iter<imx){
          n2LLold <- n2LLnew
          for(j in 1:nre){
            idx <- (1+cdm[j]):cdm[j+1]
            trval <- sum(diag(Dmat[idx,idx]))
            gval <- (rdm[j]/tau[j])-(trval/(tau[j]^2))
            gvec[j] <- (-1)*gval+crossprod(beta[idx])/((tau[j]^2)*sig)
            hmat[j,j] <- 2*(trval/(tau[j]^3)) - (rdm[j]/(tau[j]^2)) - sum(Dmat[idx,idx]^2)/(tau[j]^4)
            hmat[j,j] <- hmat[j,j] + 2*crossprod(beta[idx])/((tau[j]^3)*sig) - 2*crossprod(beta[idx],Dmat[idx,idx]%*%beta[idx])/((tau[j]^4)*sig)
            if(j>1){
              for(k in 1:(j-1)){
                kidx <- (1+cdm[k]):cdm[k+1]
                hmat[j,k] <- hmat[k,j] <- (-1)*( sum(Dmat[idx,kidx]^2) + 2*crossprod(beta[idx],Dmat[idx,kidx]%*%beta[kidx])/sig ) / ((tau[j]^2)*(tau[k]^2))
              }
            }
          }
        } else {
          ival <- FALSE
        }
        
      } # end while(ival)
      
    } else {
      # Expectation Maximization
      
      # iterative update
      while(ival) {
        
        # update tau parameter estimates
        for(j in 1:nre){
          jidx <- (1+cdm[j]):cdm[j+1]
          tau[j] <- (mean(beta[jidx]^2)/sig) + sum(diag(Dmat[jidx,jidx]))/rdm[j]
        }
        if(any(tau<=0)) tau[tau<=0] <- 10^-3
        
        # update Dmat, sig, and n2LL
        xeig <- eigen(ZtZX+diag(rep(1/tau,times=rdm)),symmetric=TRUE)
        nze <- sum(xeig$val>xeig$val[1]*.Machine$double.eps)
        Dmat <- xeig$vec[,1:nze]%*%((xeig$val[1:nze]^-1)*t(xeig$vec[,1:nze]))
        Bmat <- XtXiZ%*%Dmat
        alpha <- (XtXi+Bmat%*%t(XtXiZ))%*%Xty-Bmat%*%Zty
        beta <- Dmat%*%ZtyX
        sig <- (yty-t(c(Xty,Zty))%*%c(alpha,beta))/ndf
        if(sig<=0) sig <- 10^-3
        n2LLnew <- sum(log(xeig$val[1:nze])) + sum(log(tau)*rdm) + ndf*log(sig)
        
        # check for convergence
        vtol <- abs((n2LLold-n2LLnew)/n2LLold)
        iter <- iter + 1L
        if(vtol>tol && iter<imx){
          n2LLold <- n2LLnew
        } else {
          ival <- FALSE
        }
        
      } # end while(ival)
      
    } # end if(alg=="FS")
    
    list(tau=as.numeric(tau),sig=as.numeric(sig),iter=iter,
         cnvg=as.logical(ifelse(vtol>tol,FALSE,TRUE)),vtol=as.numeric(vtol),
         alpha=alpha,beta=beta)
    
  }