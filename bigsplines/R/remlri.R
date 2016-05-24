remlri <-
  function(yty,Xty,Zty,XtX,ZtZ,XtZ,ndf,tau=1,imx=100,tol=10^-5,alg=c("FS","NR","EM")){
    ###### REML Estimation of Random Intercept (via Fisher Scoring or EM)
    ###### Nathaniel E. Helwig (helwig@umn.edu)
    ###### Last modified: November 22, 2015
    
    ### Inputs:
    # yty: crossprod(y)
    # Xty: crossprod(X,y)
    # Zty: crossprod(Z,y)
    # XtX: crossprod(X)
    # ZtZ: diag(crossprod(Z))
    # XtZ: crossprod(X,Z)
    # ndf: nrow(X) - ncol(X)
    # tau: initial estimate of variance component
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
    #   beta is q by 1 random effect coefficients
    #     note: beta ~ N(0,sig*tau*diag(q))
    #   e is n by 1 residual vector
    #     note: e ~ N(0,sig*diag(n))
    
    ### get initial info
    XtXi <- pinvsm(XtX)
    XtXiZ <- XtXi%*%XtZ
    ZtZX <- (-1)*crossprod(XtZ,XtXiZ)
    diag(ZtZX) <- diag(ZtZX) + ZtZ
    ZtyX <- Zty-crossprod(XtZ,XtXi)%*%Xty
    nz <- length(Zty)
    Deig <- eigen(ZtZX,symmetric=TRUE)
    nze <- sum(Deig$val>Deig$val[1]*.Machine$double.eps)
    Deig$values <- Deig$values[1:nze]
    Deig$vectors <- Deig$vectors[,1:nze]
    
    ### initialize Dmat, sig, and n2LL
    newval <- 1/(Deig$val+1/tau)
    Dmat <- Deig$vec%*%tcrossprod(diag(newval),Deig$vec)
    Bmat <- XtXiZ%*%Dmat
    alpha <- (XtXi+Bmat%*%t(XtXiZ))%*%Xty - Bmat%*%Zty
    beta <- Dmat%*%ZtyX
    sig <- (yty-t(c(Xty,Zty))%*%c(alpha,beta))/ndf
    if(sig<=0) sig <- 10^-3
    n2LLold <- sum(-log(newval)) + sum(log(tau)*nz) + ndf*log(sig)
    
    ### miscellanous initializations
    alg <- alg[1]
    vtol <- 1
    iter <- 0
    ival <- TRUE
    
    ### which algorithm
    if(alg=="FS"){
      # Fisher scoring
      
      # initialize score and information
      trval <- sum(newval)
      gg <- (crossprod(beta)/((tau^2)*sig)) - ((nz/tau)-(trval/(tau^2)))
      hh <- (nz/(tau^2)) - 2*(trval/(tau^3)) + sum(newval^2)/(tau^4)
      
      # iterative update
      while(ival) {
        
        # update tau parameter estimates
        tau <- tau + gg/hh
        if(tau<=0) tau <- 10^-3
        
        # update Dmat, sig, and n2LL
        newval <- 1/(Deig$val+1/tau)
        Dmat <- Deig$vec%*%tcrossprod(diag(newval),Deig$vec)
        Bmat <- XtXiZ%*%Dmat
        alpha <- (XtXi+Bmat%*%t(XtXiZ))%*%Xty - Bmat%*%Zty
        beta <- Dmat%*%ZtyX
        sig <- (yty-t(c(Xty,Zty))%*%c(alpha,beta))/ndf
        if(sig<=0) sig <- 10^-3
        n2LLnew <- sum(-log(newval)) + sum(log(tau)*nz) + ndf*log(sig)
        
        # check for convergence (and update score and information)
        vtol <- abs((n2LLold-n2LLnew)/n2LLold)
        iter <- iter + 1L
        if(vtol>tol && iter<imx){
          n2LLold <- n2LLnew
          trval <- sum(newval)
          gg <- (crossprod(beta)/((tau^2)*sig)) - ((nz/tau)-(trval/(tau^2)))
          hh <- (nz/(tau^2)) - 2*(trval/(tau^3)) + sum(newval^2)/(tau^4)
        } else {
          ival <- FALSE
        }
        
      } # end while(ival)
      
    } else if(alg=="NR"){
      # Newton-Raphson
      
      # initialize score and information
      trval <- sum(newval)
      gg <- (crossprod(beta)/((tau^2)*sig)) - ((nz/tau)-(trval/(tau^2)))
      hh <- 2*(trval/(tau^3)) - (nz/(tau^2)) + sum(newval^2)/(tau^4)
      hh <- hh + 2*crossprod(beta)/((tau^3)*sig) - 2*crossprod(beta,Dmat%*%beta)/((tau^4)*sig)
      
      # iterative update
      while(ival) {
        
        # update tau parameter estimates
        tau <- tau + gg/hh
        if(tau<=0) tau <- 10^-3
        
        # update Dmat, sig, and n2LL
        newval <- 1/(Deig$val+1/tau)
        Dmat <- Deig$vec%*%tcrossprod(diag(newval),Deig$vec)
        Bmat <- XtXiZ%*%Dmat
        alpha <- (XtXi+Bmat%*%t(XtXiZ))%*%Xty - Bmat%*%Zty
        beta <- Dmat%*%ZtyX
        sig <- (yty-t(c(Xty,Zty))%*%c(alpha,beta))/ndf
        if(sig<=0) sig <- 10^-3
        n2LLnew <- sum(-log(newval)) + sum(log(tau)*nz) + ndf*log(sig)
        
        # check for convergence (and update score and information)
        vtol <- abs((n2LLold-n2LLnew)/n2LLold)
        iter <- iter + 1L
        if(vtol>tol && iter<imx){
          n2LLold <- n2LLnew
          trval <- sum(newval)
          gg <- (crossprod(beta)/((tau^2)*sig)) - ((nz/tau)-(trval/(tau^2)))
          hh <- 2*(trval/(tau^3)) - (nz/(tau^2)) + sum(newval^2)/(tau^4)
          hh <- hh + 2*crossprod(beta)/((tau^3)*sig) - 2*crossprod(beta,Dmat%*%beta)/((tau^4)*sig)
        } else {
          ival <- FALSE
        }
        
      } # end while(ival)
      
    } else {
      # Expectation Maximization
      
      # iterative update
      while(ival) {
        
        # update tau parameter estimates
        tau <- (mean(beta^2)/sig) + sum(newval)/nz
        if(tau<=0) tau <- 10^-3
        
        # update Dmat, sig, and n2LL
        newval <- 1/(Deig$val+1/tau)
        Dmat <- Deig$vec%*%tcrossprod(diag(newval),Deig$vec)
        Bmat <- XtXiZ%*%Dmat
        alpha <- (XtXi+Bmat%*%t(XtXiZ))%*%Xty - Bmat%*%Zty
        beta <- Dmat%*%ZtyX
        sig <- (yty-t(c(Xty,Zty))%*%c(alpha,beta))/ndf
        if(sig<=0) sig <- 10^-3
        n2LLnew <- sum(-log(newval)) + sum(log(tau)*nz) + ndf*log(sig)
        
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
         alpha=alpha,beta=beta,logLik=(-0.5)*n2LLnew)
    
  }