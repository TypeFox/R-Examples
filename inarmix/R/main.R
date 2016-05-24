inarmix <- function(formula,nclasses=1,id,data,initvals=NULL,maxiter=200,
                    stoptol=1e-5,num.restarts=1,time=NULL,dropthresh=.01)  {
  
  mod <- match.call(expand.dots=FALSE)
  modterms <- terms(formula)
  
  #### Sort the data: first by id, then by time
  subj.col <- which(colnames(data) == mod$id)
  if(!is.null(mod$time))  {
    time.col <- which(colnames(data) == mod$time)
    if(length(time.col)==0) {
      stop("time variable not found")
    }    
  }
  #### check that id variable is in the data 
  if(length(subj.col)==0) {
    stop("id variable not found")
  }    
  #### check for missing data
  if(any(is.na(data))){
    stop("There are missing data points")
  }
  if(is.null(mod$time))  {
    data <- data[order(data[,subj.col]),]
  }
  else {
    data <- data[order(data[,subj.col],data[,time.col]),]
  }
  #### Extract the responses and the model matrix
  modframe <- model.frame(formula,data)
  yy <- model.response(modframe)
  XX <- model.matrix(formula,modframe)
  
  QR <- qr(XX)
  if (QR$rank < ncol(XX))  { 
    stop("The design matrix is rank-deficient")
  }
  
  if(any(is.infinite(yy)) | any(is.infinite(XX))) {
    
  }
  #### Get num of observations per subject
  subject.id <- data[,subj.col]
  m <- length(unique(subject.id))
  nvec <- as.vector(table(subject.id))
  
  loglik.record <- numeric(num.restarts)
  initvals.record <- final.record <- list()
  
  #### Create two lists: ylist, Xlist
  #### ith component of ylist has the observations for the ith subject
  #### ith component of Xlist has the design matrix for the ith subject
  #### These lists are useful for updating the class probabilities
  ylist <- Xlist <- list()
  
  subject.names <- rep(" ",m)
  ncum <- 1 + c(0,cumsum(nvec))
  for(i in 1:m)  {
    ind1 <- ncum[i]
    ind2 <- ncum[i+1] - 1
    ylist[[i]] <- yy[ind1:ind2]
    Xlist[[i]] <- as.matrix(XX[ind1:ind2,])
    
    tmp <- unique(data[ind1,subj.col])
    subject.names[i] <- as.character(tmp)
  }
  ### or ylist <- split(yy,factor(subject.names))
  ### There is not an efficient split.matrix method 
  
  if(nclasses==1)     {
    oneclass.ests <- OneclassEsts(yy,XX,nvec,m)
    
    covs <- OuterProd(beta=oneclass.ests$beta,alpha=oneclass.ests$alpha,
                      gamma=oneclass.ests$gamma,nvec=nvec,ylist=ylist,
                      Xlist=Xlist)  
    cov.mat <- tcrossprod(solve(covs$D,covs$V),solve(covs$D))
    std.errs <- sqrt(diag(cov.mat))
    names(std.errs) <- c(colnames(XX),"autocorr.","scale")
    
    betalist <- list()
    betalist[[1]] <- oneclass.ests$beta
    tmp.pp <- PostProbs(betalist,oneclass.ests$alpha,oneclass.ests$gamma,
                        1,ylist,Xlist,m,1) 
    post.probs <- tmp.pp[[1]]
    loglikhood <- tmp.pp[[2]]
    
    results <- list() 
    class(results) <- "inarmix"
    
    results$call <- mod
    results$coefficients <- oneclass.ests$beta
    results$alpha <- oneclass.ests$alpha
    results$gamma <- oneclass.ests$gamma
    results$mix.prop <- NULL
    results$coefnames <- colnames(XX)
    results$nclasses <- 1
    results$post.probs <- NULL
    results$niter <- NULL
    results$GEE.conv <- NULL  
    results$loglikhood <- loglikhood
    results$em.converged <- NULL
    results$pss <- NULL
    results$initvals <- NULL
    results$cov.mat <- cov.mat
    results$std.errs <- std.errs
    n.pars <- length(beta) + 2 
    results$bic <- (-2)*loglikhood + n.pars*log(sum(nvec))
    results$aic <- (-2)*loglikhood + 2*n.pars
    return(results) 
  }
  else  {
    
    #### set up Matrix storage
    N <- sum(nvec)
    A1.st <- Diagonal(N,x=rep(1,N))
    A2.st <- Diagonal(N,x=rep(1,N))
    
    #### set up R.alpha "basis" matrices
    tmp <- BuildBasis(nvec)
    E1.st <- tmp$E1
    E2.st <- tmp$E2
    E3.st <- tmp$E3
    E4.st <- tmp$E4
    R.alpha <- E1.st + E2.st + E3.st + E4.st 
    
    for(rr in 1:num.restarts)  {
      #### Initialize parameter values
      if(is.null(initvals) & rr==1)  {
        initvals <- InitializePars(formula,data,subj.col,nclasses,nvec,yy,XX,ylist,Xlist)  
        
        beta <- initvals$beta
        alpha <- initvals$alpha
        gamma <- initvals$gamma
        mix.prop <- initvals$mix.prop
        
        initvals.record[[1]] = initvals
      }
      else if (rr==1 & !is.null(initvals)) {
        beta <- initvals$coef
        alpha <- initvals$autocorr
        gamma <- initvals$scale - 1
        mix.prop <- initvals$mix.prop
        
        check.inits <- (class(beta)=="matrix") & (class(alpha)=="numeric") & 
          (class(gamma)=="numeric") & (class(mix.prop)=="numeric")  
        if(!check.inits)
        {
          stop("Initial values are not in the correct format")
        }
        
        check.dims <- (nrow(beta)==nclasses) & (length(alpha)==nclasses) & 
          (length(gamma)== nclasses) & (length(mix.prop)==nclasses)
        if(!check.dims) {
          stop("Initial values do not match the number of classes")
        }
        initvals.record[[1]] <- initvals  
      }
      else if (rr > 1) {
        invals <- MultiStart(parlist,nclasses,ylist,Xlist)
        beta <- invals$beta
        alpha <- invals$alpha
        gamma <- invals$gamma
        mix.prop <- invals$mix.prop
        
        initvals.record[[rr]] <- invals
      }
      
      newbeta <- matrix(0,nrow=nrow(beta),ncol=ncol(beta))
      newalpha <- rep(0,length(alpha))
      newgamma <- rep(0,length(gamma))
      
      count <- 0
      em.converged <- 0
      current.max <- 1
      GEE.conv <- matrix(0,nrow=maxiter+1,ncol=nclasses)
      loglikhood <- pss <- rep(0,maxiter+1)
      
      done <- FALSE       
      betalist <- lapply(1:nrow(beta), function(i) beta[i,])
      tmp.pp <- PostProbs(betalist,alpha,gamma,mix.prop,ylist,Xlist,m,nclasses)
      
      post.probs <- tmp.pp$postprobs
      loglikhood[1] <- tmp.pp$loglik
      
      for (i in 1:maxiter)  {
        
        if(maxiter==0) {
          count <- 1
          break
        }
        
        #### Update parameters
        for (j in 1:nclasses)  {
          initpars <- list(beta=beta[j,],alpha=alpha[j],gamma=gamma[j])
          tmp <- GEEests(post.probs[j,],initpars,XX,yy,nvec,R.alpha,A1.st,A2.st,E1.st,E2.st,E3.st,E4.st)
          
          GEE.conv[i,j] <- tmp$conv
          
          newparams <- EnforceConstraints(tmp,Xlist,XX,subject.id)
          newbeta[j,] <- newparams$beta
          newalpha[j] <- newparams$alpha
          newgamma[j] <- newparams$gamma 
        } 
        newmix.prop <- rowSums(post.probs)/m
        
        betalist <- lapply(1:nrow(newbeta), function(i) newbeta[i,])
        tmp.pp <- PostProbs(betalist,newalpha,newgamma,newmix.prop,ylist,Xlist,m,nclasses)
        post.probs <- tmp.pp$postprobs
        loglikhood[i+1] <- tmp.pp$loglik
        
        par.stack <- c(c(newbeta),newalpha,newgamma,newmix.prop[-nclasses])
        pss.vec <- PsiFn(par.stack,len.beta=ncol(newbeta),nclasses,yy,XX,ylist,Xlist,nvec,E1.st,E2.st,E3.st,E4.st) 
        pss[i] <- sum(pss.vec^2)
        
        #### check convergence
        if(i > 1){
            done <- all(abs(pss.vec) < stoptol)
        }
        if(is.na(done)) {
          print("Warning: some post. probs are NA")
          done <- FALSE
        }
        if (done)  {
          em.converged <- 1
          count <- count + 1
          beta <- newbeta
          gamma <- newgamma
          alpha <- newalpha
          mix.prop <- newmix.prop
          break
        }
        beta <- newbeta
        gamma <- newgamma
        alpha <- newalpha
        mix.prop <- newmix.prop
        
        ### Eliminate very small classes
        if(any(mix.prop < dropthresh)) {
          if(nclasses==2) {
            print("Try a one class model")
            ## Edit later
          }
          else {
            print("Very small category, reducing the number of classes")
            ## Just remove the smallest category
            remove <- which(mix.prop==min(mix.prop))
            
            mix.prop <- mix.prop[-remove]/sum(mix.prop[-remove])
            nclasses <- nclasses - 1
            beta <- as.matrix(beta[-remove,])
            
            alpha <- alpha[-remove]
            gamma <- gamma[-remove]
            
            newbeta <- matrix(0,nrow=nclasses,ncol=ncol(beta))
            newalpha <- rep(0,nclasses)
            newgamma <- rep(0,nclasses)
          } 
        }
        cat("Iteration: ",count + 1,"\n")
        count <- count + 1
      }
      if(rr==1) { 
        numiters <- count
        EMconv <- em.converged
        Loglik <- loglikhood[1:count]
        GEEmat <- GEE.conv[1:count,]
      }
      
      parlist <- list()
      parlist$coef <- beta
      parlist$auto.corr <- alpha
      parlist$scale <- gamma + 1
      parlist$mix.prop <- mix.prop
      
      
      loglik.record[rr] <- loglikhood[count]
      final.record[[rr]] <- parlist
      
      if(rr > 1)  {
        if(loglik.record[rr] > loglik.record[rr-1]) {
          current.max <- rr
          numiters <- count
          EMconv <- em.converged
          Loglik <- loglikhood[1:count]
          GEEmat <- GEE.conv[1:count,]
        }                   
      }
      
      
      if(rr==num.restarts){
        ### pick parameters from the best run
        beta <- final.record[[current.max]]$coef
        alpha <- final.record[[current.max]]$auto.corr
        gamma <- final.record[[current.max]]$scale - 1
        mix.prop <- final.record[[current.max]]$mix.prop
        
        #### re-order according to the values of the mixture probabilities
        reord <- order(-mix.prop)
        beta <- as.matrix(beta[reord,])
        alpha <- alpha[reord]
        gamma <- gamma[reord]
        mix.prop <- mix.prop[reord]
        
        ### Update posterior probabilities using the final parameter estimates 
        betalist <- lapply(1:nrow(beta), function(i) beta[i,])
        tmp.pp <- PostProbs(betalist,alpha,gamma,mix.prop,ylist,Xlist,m,nclasses) 
        post.probs <- tmp.pp$postprobs
        
        
        q <- ncol(beta) + 2
        par.stack <- numeric(nclasses*(q + 1) - 1)
        for(i in 0:(nclasses-1))  {
          par.stack[(i*q+1):((i+1)*q)] <- c(beta[i+1,],alpha[i+1],gamma[i+1])
        }
        par.stack[(q*nclasses + 1):(nclasses*(q+1) - 1)] <- mix.prop[-nclasses]
        
        
        #### Compute Standard Errors
        len.beta <- ncol(beta)
        hh <- ComputeHess(par.stack,post.probs,len.beta,ylist,Xlist,nvec)
        if(rcond(hh$dermat) < sqrt(.Machine$double.eps)) {
          print(rcond(hh$dermat))
          hh$dermat <- hh$dermat + diag(rep(2e-8,nrow(hh$dermat)))
          print(rcond(hh$dermat))
        }
        cov.mat <- tcrossprod(solve(hh$dermat,hh$sqmat),solve(hh$dermat))
        std.errs <- sqrt(diag(cov.mat))
        
        tmp <- matrix("",ncol(XX)+2,nclasses)
        mix.names <- rep(c(""),nclasses)
        for(k in 1:nclasses)  {
          for(j in 1:ncol(XX))  {
            tmp[j,k] <- paste(colnames(XX)[j],k,sep="")
            
          }
          tmp[ncol(XX)+1,k] <- paste("autocorr.",k,sep="")
          tmp[ncol(XX)+2,k] <- paste("scale",k,sep="") 
          mix.names[k] <- paste("mix",k)
        }
        names(std.errs) <- c(c(tmp),mix.names[-nclasses])
        
        ### Transpose posterior probability matrix
        post.probs <- t(post.probs)
        rownames(post.probs) <- subject.names
        colnames(post.probs) <- rownames(beta) <- rep(" ",nclasses)
        colnames(beta) <- colnames(XX)
        for(k in 1:nclasses)   {
          colnames(post.probs)[k] <- paste("group ",k,sep="")
          rownames(beta)[k] <- paste("group ",k,sep="")
        }
        
        results <- list() 
        class(results) <- "inarmix"
        
        results$call <- mod
        results$coefficients <- beta
        results$alpha <- alpha
        results$gamma <- gamma
        results$mix.prop <- mix.prop
        results$coefnames <- colnames(XX)
        results$post.probs <- post.probs
        results$fitted.values <- exp(tcrossprod(XX,beta))
        res_se <- sqrt(t(t(results$fitted.values)*(gamma + 1)))
        results$residuals <- (yy - results$fitted.values)/res_se
        results$niter <- numiters
        results$GEE.conv <- GEEmat  
        results$loglikhood <- Loglik
        results$em.converged <- EMconv
        results$pss <- pss[1:count]
        results$initvals <- initvals.record[[current.max]]
        results$nclasses <- nclasses
        results$cov.mat <- cov.mat
        results$std.errs <- std.errs
        n.pars <- nclasses*(ncol(beta) + 3) - 1
        results$bic <- (-2)*Loglik[numiters] + n.pars*log(sum(nvec))
        results$aic <- (-2)*Loglik[numiters] + 2*n.pars
        
        if(num.restarts==1)  {
          results$startingvals <- NULL
          results$reploglik <- NULL
          results$finalvals <- NULL
        }
        else {
          results$startingvals <- initvals.record
          results$reploglik <- loglik.record
          results$finalvals <- final.record
        }
      }
      else {
        cat("\n Replication:",rr + 1,"\n")
      }
    }
    return(results)    
  }
}

## Call to C++ for the posterior probabilities
PostProbs <- function(betalist,alpha,gamma,mix,ylist,Xlist,m,nclasses) {
  .Call("PostProb",betalist,alpha,gamma,mix,ylist,Xlist,m,nclasses,PACKAGE="inarmix")
}

## Call to C++ for computing the estimated covariance matrix
ComputeHess <- function(parstack,postprobs,lenbeta,ylist,Xlist,nvecs) {
  .Call("ComputeHessian",parstack, postprobs, lenbeta, ylist, Xlist, nvecs, PACKAGE="inarmix")
}
