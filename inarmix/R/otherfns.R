predict.class <- function(results)  {
  ### results should be an object of class inarmix
  a <- apply(results$post.probs,1,function(x) which(x==max(x))[1])
  return(a)
}

EnforceConstraints <- function(params,Xlist,XX,subj.factor)  {    
  beta <- as.vector(params$beta)
  ff <- function(X) { return(X %*% beta)}
  gg <- function(x) { 
    tmp <- abs(diff(x))
    ans <- ifelse(length(tmp) > 0,max(tmp),0)
    return(ans)
  }
  lin.stack <- XX%*%beta
  
  split.lin <- split(lin.stack,subj.factor)
  bounds <- unlist(lapply(split.lin,gg))
  
  bound <- max(bounds)
  alpha.bound <- exp(-bound/2) + 1e-4
  newalpha <- min(params$alpha,alpha.bound)
  
  alpha <- max(newalpha,1e-4)
  alpha <- min(alpha,1 - 1e-4)
  gamma <- max(params$gamma,1e-4)
  
  return(list(beta=beta,alpha=alpha,gamma=gamma))
}


MyPoisson_Reg <- function(yy,XX,numiter=5)  {
  ### This function computes regression coefficients for
  ### a Poisson regression under a log-link function.
  
  d <- ncol(XX)
  beta <- rep(1,d) 
  
  ### Perform 5 iterations of the Fisher 
  ### scoring algorithm
  for(k in 1:numiter)  {
    mu <- as.vector(exp(XX%*%beta))
    W <- diag(mu)
    A <- crossprod(W%*%XX,XX)
    
    ### Check that the condition number of the
    ### above matrix is not extremely large
    if(rcond(A) < sqrt(.Machine$double.eps)) {
      A <- A + diag(rep(2e-8,d))
    }
    ### update beta
    beta <- beta + solve(A,crossprod(XX,yy-mu))
  }
  return(beta)
}





InitializePars = function(formula,data,subj.col,nclasses,nvec,yy,XX,ylist,Xlist)  {
  
  beta <- matrix(0,nrow=nclasses,ncol=ncol(XX))
  alpha <- gamma <- mix.prop <- rep(1/nclasses,nclasses)
  deviance <- matrix(0,nrow=length(yy),ncol=nclasses)
  m <- length(nvec)
  fitted <- numeric(length(yy))
  subject.ind <- rep(1:m,nvec)
  
  #### Sample nclasses number of people
  subject.init <- sample(1:m,size=nclasses)
  for(k in 1:nclasses)  {
    #tmpdata <- data.frame(y=ylist[[subject.init[k]]],Xlist[[subject.init[k]]])
    current_y <- ylist[[subject.init[k]]]
    current_X <- Xlist[[subject.init[k]]]
    ### run glm 
    ### still need to worry if there are fewer observations
    ### than predictors.
    if(length(current_y) < ncol(XX)) {
      ss <- c(1:m)
      orig.subj <- subject.init[[k]]
      add.subj <- sample(ss[-orig.subj],size=1)
      ty <- unlist(ylist[c(orig.subj,add.subj)])
      tx <- rbind(Xlist[[orig.subj]],Xlist[[add.subj]])
      
      #tmpdata <- data.frame(y=ty,tx)
      if(length(ty) >= ncol(XX)) {
        #suppressWarnings(bet <- glm(y ~ . - 1,data=tmpdata,family=poisson(),maxit=5)$coef)
        bet <- MyPoisson_Reg(ty,tx)
      }
      else {
        bet <- rep(1,ncol(XX))
      }
    }
    else {
      # suppressWarnings(bet <- glm(y ~ . - 1,data=tmpdata,family=poisson(),maxit=5)$coef)
      bet <- MyPoisson_Reg(current_y,current_X)             
    }
    fitted <- exp(XX%*%bet)
    deviance[,k] <- 2*fitted - 2*yy*log(fitted)
  }
  tau <- 3
  for(i in 1:tau)   {
    ### make assignments based on minimum deviance.
    RR <- aggregate(deviance,by=list(subject.ind),FUN="sum")
    assign <- apply(RR[,-1],1,function(x) which.min(x))
    assign.expand <- rep(assign,nvec)
    
    tt <- data.frame(subject=subject.ind,class=assign.expand)
    for(k in 1:nclasses)  {
      
      assign.ind <- which(tt$class==k)
      
      if(length(assign.ind)==0)  {
        assign.ind <- sample(1:m,size=floor(m/nclasses))
      }              
      
      tmp.yy <- as.vector(yy[assign.ind])
      tmp.XX <- as.matrix(XX[assign.ind,])
      
      if(i==tau)  {
        tmp.nvec <- as.vector(nvec[assign==k])
        tmp.m <- length(tmp.nvec)
        
        myfit <- OneclassEsts(yy=tmp.yy,XX=tmp.XX,nvec=tmp.nvec,m=tmp.m)
        beta[k,] <- as.numeric(myfit$beta)
        
        alpha[k] <- max(myfit$alpha,.05)
        gamma[k] <- max(myfit$gamma,.1)
      }
      else {
        if(length(assign.ind) == 0)  {
          bet <- bet + rnorm(length(bet))
        }
        else {
          bet <- MyPoisson_Reg(tmp.yy,tmp.XX) 
          
          fitted <- exp(XX%*%bet)
          deviance[,k] <- 2*fitted - 2*yy*log(fitted)
        }
      }
    }
  }
  ### In case any NA's are present
  beta[which(is.na(beta))] <- 1
  alpha[which(is.na(alpha))] <- .2
  gamma[which(is.na(gamma))] <- 1
  
  tmp.mix <- table(assign)/m
  if(length(tmp.mix) < nclasses) {
    a <- rep(1,nclasses-length(tmp.mix))
    mix.prop <- c(tmp.mix,a)/(sum(tmp.mix) + sum(a))
  }
  else {
    mix.prop <- tmp.mix
  }
  return(list(beta=beta,alpha=alpha,gamma=gamma,mix.prop=mix.prop))
}


MultiStart <- function(parlist,nclasses,ylist,Xlist)  {
  
  m <- length(ylist)
  alpha <- gamma <- mix <- numeric(nclasses)
  beta.mat <- matrix(0,nrow=nclasses,ncol=ncol(Xlist[[1]]))
  
  aa <- rnorm(nclasses)
  subject.ind <- sample(1:m,size=nclasses)
  for(k in 1:nclasses)  {
    #tmpdata <- data.frame(y=ylist[[subject.ind[k]]],Xlist[[subject.ind[k]]])
    #suppressWarnings(beta.vals <- glm(y ~ . - 1,data=tmpdata,family="poisson"))
    
    beta.vals <- MyPoisson_Reg(yy=ylist[[subject.ind[k]]],XX=Xlist[[subject.ind[k]]]) 
    resids <- ylist[[subject.ind[k]]] - exp(Xlist[[subject.ind[k]]]%*%beta.vals)
    
    beta.mat[k,] <- beta.vals  
    gamma[k] <- max(mean(resids^2) - 1,.1)
    alpha[k] <- runif(1)
    mix[k] <- exp(aa[k])/sum(exp(aa))
  }
  
  ### Take a convex combination of previous values and new random values
  ### for the parameters
  beta.mat <- .75*parlist$coef + .25*beta.mat
  alpha <- .75*parlist$auto.corr + .25*alpha
  gamma <- .75*(parlist$scale - 1) + .25*gamma
  mix <- .75*parlist$mix.prop + .25*mix
  
  inpars <- list()
  inpars$beta <- beta.mat
  inpars$alpha <- alpha
  inpars$gamma <- gamma
  inpars$mix.prop <- mix
  
  return(inpars)
}


PsiFn <- function(par.stack,len.beta,nclasses,yy,XX,ylist,Xlist,nvec,E1,E2,E3,E4,group.remove=NULL)  {
  
  beta.vec <- par.stack[1:(len.beta*nclasses)]
  beta <- matrix(beta.vec,nrow=nclasses,ncol=len.beta)
  alpha.start <- len.beta*nclasses + 1 
  alpha <- par.stack[alpha.start:(alpha.start + nclasses - 1)]
  gam.start <- alpha.start + nclasses
  gamma <- par.stack[gam.start:(gam.start + nclasses - 1)]
  mix.start <- gam.start + nclasses
  tmp.mix <- par.stack[mix.start:(mix.start + nclasses - 2)]
  mix.prop <- c(tmp.mix,1 - sum(tmp.mix))
  
  est.eq <- rep(0,length(par.stack))
  
  m <- length(nvec)
  betalist <- lapply(1:nrow(beta), function(i) beta[i,])
  tmp.pp <- PostProbs(betalist,alpha,gamma,mix.prop,ylist,Xlist,m,nclasses) 
  post.probs <- tmp.pp$postprobs
  
  N <- sum(nvec)
  m <- length(nvec)
  A1.store <- Diagonal(N,x=rep(1,N))
  A2.store <- Diagonal(N,x=rep(1,N))
  ##### initialize R.alpha, and A.store
  
  ind1 <- 0
  ind2 <- 0
  for(k in 1:nclasses)   {
    ind1 <- ind2 + 1
    ind2 <- ind2 + len.beta
    weights <- post.probs[k,]
    expanded.weights <- rep(weights,nvec)
    
    mu.stack <- as.vector(exp(XX%*%beta[k,]))
    diag(A1.store) <- sqrt(expanded.weights*mu.stack)
    diag(A2.store) <- (1/sqrt(mu.stack))
    R.alpha <- (alpha[k]*E1 + (1 + alpha[k]^2)*E2 + E3)/(1 - alpha[k]^2) + E4 
    weighted.dat <- sqrt(expanded.weights)*(yy - mu.stack)    
    est.eq[ind1:ind2] <- as.vector(crossprod(A1.store%*%XX,R.alpha%*%A2.store%*%weighted.dat))
    
  }
  alpha.ind <- len.beta*nclasses + 1
  gamma.ind <- alpha.ind + nclasses
  mix.ind <- seq(gamma.ind + nclasses,gamma.ind + 2*nclasses-2,by=1)
  for(k in 1:nclasses)  {
    weights <- post.probs[k,]
    expanded.weights <- rep(weights,nvec)
    mu.stack <- as.vector(exp(XX%*%beta[k,]))
    diag(A1.store) <- (1/sqrt(mu.stack))
    weighted.dat <- sqrt(expanded.weights)*(yy - mu.stack)
    
    #### alpha est.eq
    #### First update R.alpha as R.sandwich
    N.weight <- sum(weights*(nvec - 1))
    R.alpha <-  -((1 + alpha[k]^2)*E1 + 4*alpha[k]*E2 + 2*alpha[k]*E3)/((1 - alpha[k]^2)^2)
    LHT <- as.numeric(crossprod(A1.store%*%weighted.dat,R.alpha%*%A1.store%*%weighted.dat))
    RHT <- (2*alpha[k]*N.weight*(1+gamma[k]))/(1 - alpha[k]^2)
    est.eq[alpha.ind] <- LHT + RHT
    
    #### gamma est.eq
    #### update R.alpha
    N.weight <- sum(weights*nvec)
    R.alpha <- (alpha[k]*E1 + (1 + alpha[k]^2)*E2 + E3)/(1 - alpha[k]^2) + E4 
    
    comp2 <- (1/(1 + gamma[k]))*as.numeric(crossprod(A1.store%*%weighted.dat,R.alpha%*%A1.store%*%weighted.dat)) 
    est.eq[gamma.ind] <- comp2 - N.weight
    
    gamma.ind <- gamma.ind + 1
    alpha.ind <- alpha.ind + 1
  }
  est.eq[mix.ind] <- rowSums(post.probs)[1:(nclasses-1)] - m*tmp.mix
  
  return(est.eq/m)
}