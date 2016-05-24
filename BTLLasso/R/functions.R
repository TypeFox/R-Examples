

loglik <- function(coef,y,design,kat){
  eta <- matrix(design%*%coef,ncol=kat-1,byrow=TRUE)
  pi.help <- matrix(exp(eta)/(1+exp(eta)),ncol=kat-1)
  # print(pi.help)
  
  pi <- pi.help
  if(kat>2){
    for (i in 2:(kat-1)){
      pi[,i] <- pi.help[,i]-pi.help[,i-1]
    }
  }
  pi <- cbind(pi,1-pi.help[,kat-1])
  
  
  yhelp <- rep(y,each=kat)
  yhelp <- matrix(as.numeric(yhelp==rep(1:kat,length(y))),byrow=T,ncol=kat)
  
  
  loglik <- sum(yhelp*log(pi))
  loglik
}


create.design <- function(X, m){
  p <- ncol(X)
  n <- nrow(X)
  
  I <- m*(m-1)/2
  
  pos1 <- 1
  pos2 <- 2
  line <- 1
  
  design.help <- matrix(0,nrow=I,ncol=m)  
  
  while(line<=nrow(design.help)){
    
    design.help[line,pos1] <- 1
    design.help[line,pos2] <- -1
    
    line <- line + 1
    pos2 <- pos2 + 1
    
    if(pos2>m){
      pos1 <- pos1 + 1
      pos2 <- pos1 + 1
    }
  }
  design.help <- design.help[,-m]
  rep.vec <- rep(1:I,each=n)
  design.help <- design.help[rep.vec,]
  
  
  design <- matrix(c(apply(cbind(1,X),2,function(xx){xx*c(design.help)})),nrow=nrow(design.help))
  
  return(design)
}




predict.BTLLasso <- function(X,coef,q,m){
  k <- q+1
  I <- choose(m,2)
  
  design2 <- create.design(X,m)
  design <- t(matrix(rep(c(design2),each=q),nrow=ncol(design2),byrow=TRUE))
  
  n.theta <- floor(q/2)
  
  if(k>2){
    theta.design <- matrix(0,ncol=n.theta,nrow=nrow(design))
    
    for(i in 1:n.theta){
      vec1 <- rep(0,q)
      vec1[c(i,q-i+1)] <- c(1,-1)
      theta.design[,i] <- rep(vec1,nrow(design2))
    }
    
    design <- cbind(theta.design, design)
  }
  
  
  eta <- matrix(design%*%coef,ncol=q,byrow=TRUE)
  pi.help <- matrix(exp(eta)/(1+exp(eta)),ncol=q)
  # print(pi.help)
  
  pi <- pi.help
  if(k>2){
    for (i in 2:(q)){
      pi[,i] <- pi.help[,i]-pi.help[,i-1]
    }
  }
  pi <- cbind(pi,1-pi.help[,q])
  
  c(t(pi))
}


reparam <- function(x){
  z <- ncol(x)+1
  K <- matrix((-1/z), ncol = z - 1, nrow = z - 1)
  diag(K) <- (z - 1)/z
  x.sym <- x %*% K
  
  
  x2 <- x.sym %*% matrix(rep(-1, z - 1), nrow= z - 1) 
  
  x.sym <- cbind(x2,x.sym)
  
  x.sym
}


expand.coefs <- function(coef, n.theta, m, symmetric = TRUE){
  
  intercepts <- coef[,(n.theta+1):(n.theta+m-1)]
  gamma <-  coef[,(n.theta+m):ncol(coef)]
  
  p <- ncol(gamma)/(m-1)
  
  if(symmetric){
    intercepts.new <- reparam(intercepts)
  }else{
    intercepts.new <- cbind(0,intercepts)
  }
  
  coef.new <- c()
  
  if(n.theta>0){
    coef.new <- coef[,1:n.theta,drop=FALSE]
  }
  
  coef.new <- cbind(coef.new,intercepts.new)
  
  
  index <- 1
  for(i in 1:p){
    if(symmetric){
      coef.new <- cbind(coef.new,reparam(gamma[,index:(index+m-2)]))
    }else{
      coef.new <- cbind(coef.new,0,gamma[,index:(index+m-2)])
    }
    index <- index+m-1
    
  }
  
  
  coef.new
}
