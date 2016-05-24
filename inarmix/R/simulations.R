GenerateMixData = function(m,coefs,autocorr,scale=NULL,mix,design.mat,return.labels=F,poisson=F)  {
  T <- nrow(design.mat)
  X <- design.mat
  
  y <- matrix(0,nrow=m,ncol=T)
  labels <- matrix(0,nrow=m,ncol=2)
  colnames(labels) <- c("subject","group label")
  nclasses <- length(mix)
  beta <- coefs
  alpha <- autocorr
  
  if(any((alpha >= 1) | (alpha <=0))) {
    stop("All autocorrelation parameters must be between zero and one")
  }
  
  if(!poisson)  {
    ### negative binomial case
    if(is.null(scale))  {
      stop("scale must be a vector of dispersion parameters")
    }
    if(any(scale <= 1))  {
      stop("All scale parameters must be greater than one")
    }
    gamma <- scale - 1
    
    for (i in 1:m)  {
      
      mixind <- sample(c(1:nclasses), size=1, prob = mix)
      labels[i,1] <- i
      labels[i,2] <- mixind
      
      theta <- exp(X%*%beta[mixind,])/gamma[mixind]
      
      theta.lag = theta[1:(T-1)];  theta.shift = theta[2:T]
      theta.c = alpha[mixind]*sqrt(theta.lag*theta.shift)
      theta.io = theta.lag - theta.c
      theta.in = theta.shift - theta.c
      
      y[i,1] = rnbinom(1,size = theta[1],prob=1/(1+gamma[mixind]))
      for (t in 2:T)  {
        I = rnbinom(1,size=theta.in[t-1],prob=1/(1+gamma[mixind]))
        u = rbeta(1,shape1 = theta.c[t-1],shape2=theta.io[t-1])
        a = rbinom(1,size=y[i,t-1],prob=u)
        y[i,t] = a + I
      }
    }
  }   
  else  {
    #### Poisson case
    
    for (i in 1:m)  {
      
      mixind <- sample(c(1:nclasses), size=1, prob = mix)
      labels[i,1] <- i
      labels[i,2] <- mixind
      
      mu <- exp(X%*%beta[mixind,])
      
      y[i,1] = rpois(1,lambda=mu[1])
      for (t in 2:T)  {
        parm1 <- mu[t] - alpha*sqrt(mu[t]*mu[t-1])
        parm2 <- alpha*sqrt(mu[t]/mu[t-1])
        
        I = rpois(1,lambda=parm1)
        W = rbinom(1,size=y[i,t-1],prob=parm2)
        y[i,t] = W + I
      }
    }
  }
  
  yy <- c(t(y))
  nc <- ncol(X)
  XX <- matrix(0,nrow=T*m,ncol=nc)
  for(k in 1:nc)  {
    XX[,k] <- rep(X[,k],m)
  }
  subj <- rep(1:m,each=T)
  res <- data.frame(cbind(yy,subj,XX))
  colnames(res) <- c("y","subject",colnames(design.mat))
  
  if(return.labels)   {
    return(list(dataset=res,labels=labels))
  }
  else  {
    return(res)
  }
}