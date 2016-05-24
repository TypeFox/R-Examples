blockstand <- function(x, ipen.which, inotpen.which)
{
  ## Author: Lukas Meier, Date:  4 Aug 2006, 08:50
  
  n <- nrow(x)
  x.ort <- x
  scale.pen <- list(); length(scale.pen) <- length(ipen.which)
  scale.notpen <- NULL
  
  if(length(inotpen.which) > 0){
    one <- rep(1, n)
    scale.notpen <- sqrt(drop(one %*% (x[,inotpen.which]^2)) / n)
    x.ort[,inotpen.which] <- scale(x[,inotpen.which], FALSE, scale.notpen)
  }
  
  for(j in 1:length(ipen.which)){
    ind <- ipen.which[[j]]
    decomp <- qr(x[,ind])
    if(decomp$rank < length(ind)) ## Warn if block has not full rank
      stop("Block belonging to columns ", paste(ind, collapse = ", "),
           " has not full rank! \n")
    scale.pen[[j]] <- qr.R(decomp) * 1 / sqrt(n)
    x.ort[,ind] <- qr.Q(decomp) * sqrt(n)
  }
  list(x = x.ort, scale.pen = scale.pen, scale.notpen = scale.notpen)
}

taylor.opt<-function(t_opt,y,X,fixef,ranef,Grad,family,P, K=NULL)
{
  delta<-c(fixef,ranef)+t_opt*Grad
  Eta<- X%*%delta
  if(is.null(K)){
    mu<-family$linkinv(Eta)
  }else{
    Eta_cat <- matrix(Eta, byrow = TRUE, ncol = K)
    mu <- c(t(family$linkinv(Eta_cat)))
  }
  
  loglik <- logLik.glmmLasso(y=y,mu=mu,family=family,ranef.logLik=NULL,penal=FALSE,K=K) 
  
  loglik<- -loglik + 0.5*t(delta[(length(fixef)+1):length(delta)]) %*% (P %*% delta[(length(fixef)+1):length(delta)])
  return(loglik)
}

taylor.opt.noRE<-function(t_opt,y,X,fixef,Grad,family, K=NULL)
{
  
  delta<-fixef+t_opt*Grad
  Eta<- X%*%delta
  if(is.null(K)){
    mu<-family$linkinv(Eta)
  }else{
    Eta_cat <- matrix(Eta, byrow = TRUE, ncol = K)
    mu <- c(t(family$linkinv(Eta_cat)))
  }
  loglik <- -logLik.glmmLasso(y=y,mu=mu,family=family,ranef.logLik=NULL,penal=FALSE, K = K) 
  
  return(loglik)
}


gradient.lasso.block<-function(score.beta,b,lambda.b,block)       
{
  p<-length(b)
  grad.lasso<-rep(0,p)
  
  block2 <- rep(block,block)
  lambda_vec<-rep(lambda.b,length(b))
  lambda_vec <- lambda_vec*sqrt(block2)
  
  group.sum<-rep(0,length(block))
  group.sum[1]<-sqrt(sum(b[1:block[1]]^2))
  if(length(block)>1)
  {  
    for (i in 2:length(block))
      group.sum[i]<-sqrt(sum(b[(sum(block[1:(i-1)])+1):sum(block[1:i])]^2))
  }
  if(group.sum[1]!=0)
  {
    grad.lasso[1:block[1]]<-score.beta[1:block[1]]-lambda_vec[1]*(b[1:block[1]]/group.sum[1])
  }else{
    
    if(group.sum[1]==0 & sqrt(sum(score.beta[1:block[1]]^2))>lambda_vec[1] )
    {
      grad.lasso[1:block[1]]<-score.beta[1:block[1]]-lambda_vec[1]*(score.beta[1:block[1]]/sqrt(sum(score.beta[1:block[1]]^2)))
    }else{
      grad.lasso[1:block[1]]<-0
    }}
  
  if(length(block)>1)
  {  
    for (i in 2:length(block))
    {
      if(group.sum[i]!=0)
      {
        grad.lasso[(sum(block[1:(i-1)])+1):sum(block[1:i])]<-score.beta[(sum(block[1:(i-1)])+1):sum(block[1:i])]-lambda_vec[sum(block[1:(i-1)])+1]*(b[(sum(block[1:(i-1)])+1):sum(block[1:i])]/group.sum[i])
      }else{
        
        if(group.sum[i]==0 & sqrt(sum(score.beta[(sum(block[1:(i-1)])+1):sum(block[1:i])]^2))>lambda_vec[sum(block[1:(i-1)])+1])
        {
          grad.lasso[(sum(block[1:(i-1)])+1):sum(block[1:i])]<-score.beta[(sum(block[1:(i-1)])+1):sum(block[1:i])]-lambda_vec[sum(block[1:(i-1)])+1]*(score.beta[(sum(block[1:(i-1)])+1):sum(block[1:i])]/sqrt(sum(score.beta[(sum(block[1:(i-1)])+1):sum(block[1:i])]^2)))
        }else{
          grad.lasso[(sum(block[1:(i-1)])+1):sum(block[1:i])]<-0
        }}
    }}
  return(grad.lasso)
}

gradient.lasso<-function(score.beta,b,lambda.b)
{
  p<-length(b)
  grad.lasso<-rep(0,p)
  
  b.isnt.0<-b!=0
  grad.lasso[b.isnt.0]<-score.beta[b.isnt.0]-lambda.b*sign(b[b.isnt.0])
  b.is.0<-(b==0 & abs(score.beta)>=lambda.b)
  grad.lasso[b.is.0]<-score.beta[b.is.0]-lambda.b*sign(score.beta[b.is.0])
  b.is.0<-(b==0 & abs(score.beta)<lambda.b)
  grad.lasso[b.is.0]<-0
  return(grad.lasso)
}


t.change<-function(grad,b)
{
  a<-(sign(b)==-sign(grad)&sign(grad)!=0)
  rate <- rep(Inf, length(a))
  rate[a] <--b[a]/grad[a]
  ret.obj<-list()
  ret.obj$min.rate<-min(rate)
  ret.obj$whichmin <- which.min(rate)
  return(ret.obj)
}

