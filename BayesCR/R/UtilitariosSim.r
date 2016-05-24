PearsonVII <- function(y,mu,sigma2,nu,delta)
{
  Acum <- z <- vector(mode = "numeric", length = length(y))
  sigma2a <- sigma2*(delta/nu)
  z <- (y-mu)/sqrt(sigma2a)
  Acum <- pt(z,df=nu)
  return(Acum)
}

AcumSlash <- function(y,mu,sigma2,nu)
{
  Acum <- z <- vector(mode = "numeric", length = length(y))
  z <- (y-mu)/sqrt(sigma2)
  for (i in 1:length(y))
  { 
    f1 <- function(u) nu*u^(nu-1)*pnorm(z[i]*sqrt(u))
    Acum[i]<- integrate(f1,0,1)$value	 	
  }
  return(Acum)
}

AcumNormalC <- function(y,mu,sigma2,nu)
{
  Acum <- vector(mode = "numeric", length = length(y))
  eta  <- nu[1]    
  gama <- nu[2]
  Acum <- eta*pnorm(y,mu,sqrt(sigma2/gama)) + (1-eta)*pnorm(y,mu,sqrt(sigma2))
  return(Acum)
}

dPearsonVII<- function(y,mu,sigma2,nu,delta)
{
  f <-  z <- vector(mode = "numeric", length = length(y))
  sigma2a <- sigma2*(delta/nu)
  z <- (y-mu)/sqrt(sigma2a)
  f <- dt(z,df=nu)/sqrt(sigma2a)
  return(f)
}

dSlash <- function(y,mu,sigma2,nu)
{
   resp <- z <- vector(mode = "numeric", length = length(y))
   z <- (y-mu)/sqrt(sigma2)
   for (i in 1:length(y)) 
   {
     f1 <- function(u) nu*u^(nu-0.5)*dnorm(z[i]*sqrt(u))/sqrt(sigma2)
     resp[i] <- integrate(f1,0,1)$value	 	
   }
     return(resp)
}

dNormalC <- function(y,mu,sigma2,nu)
{
   Acum <- vector(mode = "numeric", length = length(y))
   eta  <- nu[1]    
   gama <- nu[2]
   Acum <- eta*dnorm(y,mu,sqrt(sigma2/gama)) + (1-eta)*dnorm(y,mu,sqrt(sigma2))
   
   return(Acum)
}

verosCN <- function(auxf,cc,cens,sigma2,nu)
{
  auxf1 <- sqrt(auxf)
  if (cens=="1")
  {
    ver1 <-(sum(log(dNormalC(auxf1[cc==0],0,1,nu)/sqrt(sigma2)))+ sum(log(AcumNormalC(auxf1[cc==1],0,1,nu))))
  }
  if (cens=="2")
  {
    ver1 <-(sum(log(dNormalC(auxf1[cc==0],0,1,nu)/sqrt(sigma2)))+ sum(log(AcumNormalC(auxf1[cc==1],0,1,nu))))
  }
  return(ver1)
}

Rhat1 <- function(param,n.iter,burnin,n.chains,n.thin) 
{
  param <- as.matrix(param)
  efect <- (n.iter-burnin)/n.thin
  p <- ncol(param)
  mat <- matrix(param,nrow=efect,ncol=n.chains*p)    
  rhat <- matrix(0,nrow=p,ncol=1)                                          
  for(i in 1:p)
  {
    l1 <- 2*(i-1)+1
    c1 <- 2*i
    rhat[i,1] <- Rhat(mat[,l1:c1]) 
  }
  return(rhat=rhat)
}


Rhat <- function(mat) 
{
  m <- ncol(mat)
  n <- nrow(mat)
  b <- apply(mat,2,mean)
  B <- sum((b-mean(mat))^2)*n/(m-1)
  w <- apply(mat,2,var)
  W <- mean(w)
  s2hat <- (n-1)/n*W + B/n
  Vhat <- s2hat + B/m/n 
  covWB <- n /m * (cov(w,b^2)-2*mean(b)*cov(w,b))
  varV <- (n-1)^2 / n^2 * var(w)/m + (m+1)^2 / m^2 / n^2 * 2*B^2/(m-1) + 2 * (m-1)*(n-1)/m/n^2 * covWB
  df <- 2 * Vhat^2 / varV
  R <- sqrt((df+3) * Vhat / (df+1) / W)
  return(R)
}
 

hpd <- function(x, alpha)
{
   n <- length(x)
   m <- max(1, ceiling(alpha * n))
   y <- sort(x)
   a <- y[1:m]
   b <- y[(n - m + 1):n]
   i <- order(b - a)[1]
   structure(c(a[i], b[i]), names = c("Lower Bound", "Upper Bound"))
}

MHnu<-function(last,U,lambda,prior="Jeffreys",hyper)
{
  n <- length(U) 
  if(prior=="Jeffreys")
  {
    gJeffreys <- function(nu,U) 
    {
      n<-length(U)
      ff<- log(sqrt(nu/(nu+3))*sqrt(trigamma(nu/2)-trigamma((nu+1)/2)-2*(nu+3)/((nu)*(nu+1)^2)))+0.5*n*nu*log(nu/2)+(0.5*nu)*sum(log(U)-U)-n*log(gamma(nu/2))
      return(ff)
    }
    Fonseca1 <-deriv(~log(sqrt(nu/(nu+3))*sqrt(trigamma(nu/2)-trigamma((nu+1)/2)-2*(nu+3)/((nu)*(nu+1)^2))) +0.5*n*nu*log(nu/2)-n*log(gamma(nu/2)),c("nu"),function(nu){},hessian=TRUE)
    Fonseca2 <- deriv(~(0.5*nu)*(log(U)-U),c("nu"),function(U,nu){},hessian=TRUE)
    aux1 <- Fonseca1(last)
    aux2 <- Fonseca2(U,last)                                                   
    q1 <- attr(aux1,"gradient")[1]+sum(attr(aux2,"gradient"))
    q2 <- attr(aux1,"hessian")[1]+sum(attr(aux2,"hessian"))
    aw <- last-q1/q2
    bw <- max(0.001,-1/q2)
    cand <- rtrunc(1, spec="norm",a=2.1, b=100, mean = aw, sd = sqrt(bw)) 
    alfa <- (exp(gJeffreys(cand,U))/exp(gJeffreys(last,U)))*(dtrunc(last,spec="norm",a=2.1, b=100, mean = aw, sd = sqrt(bw))/dtrunc(cand,spec="norm", a=2.1, b=100, mean = aw, sd = sqrt(bw)))
  }
  
  if(prior=="Exp") 
  {
    gExp <- function(nu,U,hyper) 
    {
      n<-length(U)
      ff<- (-hyper*nu) + 0.5*n*nu*log(nu/2)+(0.5*nu)*sum(log(U)-U)-n*log(gamma(nu/2))
      return(ff)
    }
    Fonseca1<- deriv(~ (-hyper*nu) + 0.5*n*nu*log(nu/2)-n*log(gamma(nu/2)),c("nu"),function(nu){},hessian=TRUE)
    Fonseca2<- deriv(~(0.5*nu)*(log(U)-U),c("nu"),function(U,nu){},hessian=TRUE)
    aux1<-Fonseca1(last)
    aux2<-Fonseca2(U,last)                                                   
    q1<-attr(aux1,"gradient")[1]+sum(attr(aux2,"gradient"))
    q2<-attr(aux1,"hessian")[1]+sum(attr(aux2,"hessian"))
    aw<- last-q1/q2
    bw<- max(0.001,-1/q2)
    cand <- rtrunc(1, spec="norm",a=2.1, b=100, mean = aw, sd = sqrt(bw))                  
    alfa<-(exp(gExp(cand,U,hyper))/exp(gExp(last,U,hyper)))*(dtrunc(last,spec="norm",a=2.1, b=100, mean = aw, sd = sqrt(bw))/dtrunc(cand,spec="norm", a=2.1, b=100, mean = aw, sd = sqrt(bw)))
  }
  
  if(prior=="Unif") 
  {
    gUnif <- function(nu,U) 
    {
      n<-length(U)
      ff<- 0.5*n*nu*log(nu/2)+(0.5*nu)*sum(log(U)-U)-n*log(gamma(nu/2))
      return(ff)
    }
    Fonseca1<- deriv(~ 0.5*n*nu*log(nu/2)-n*log(gamma(nu/2)),c("nu"),function(nu){},hessian=TRUE)
    Fonseca2<- deriv(~(0.5*nu)*(log(U)-U),c("nu"),function(U,nu){},hessian=TRUE)
    aux1<-Fonseca1(last)
    aux2<-Fonseca2(U,last)                                                   
    q1<-attr(aux1,"gradient")[1]+sum(attr(aux2,"gradient"))
    q2<-attr(aux1,"hessian")[1]+sum(attr(aux2,"hessian"))
    aw<- last-q1/q2
    bw<- max(0.001,-1/q2)
    cand <- rtrunc(1, spec="norm",a=2.1, b=100, mean = aw, sd = sqrt(bw))                                  
    alfa<-(exp(gUnif(cand,U))/exp(gUnif(last,U)))*(dtrunc(last,spec="norm",a=2.1, b=100, mean = aw, sd = sqrt(bw))/dtrunc(cand,spec="norm", a=2.1, b=100, mean = aw, sd = sqrt(bw)))
  }
  
  if(prior=="Hierar") 
  {
    gHierar <- function(nu,U) 
    {
      n<-length(U)
      ff<- (-lambda*nu) + 0.5*n*nu*log(nu/2)+(0.5*nu)*sum(log(U)-U)-n*log(gamma(nu/2))
      return(ff)
    }
    Fonseca1<- deriv(~ (-lambda*nu) + 0.5*n*nu*log(nu/2)-n*log(gamma(nu/2)),c("nu"),function(nu){},hessian=TRUE)
    Fonseca2<- deriv(~(0.5*nu)*(log(U)-U),c("nu"),function(U,nu){},hessian=TRUE)
    aux1<- Fonseca1(last)
    aux2<-Fonseca2(U,last)                                                   
    q1<-attr(aux1,"gradient")[1]+sum(attr(aux2,"gradient"))
    q2<-attr(aux1,"hessian")[1]+sum(attr(aux2,"hessian"))
    aw<- last-q1/q2
    bw<- max(0.001,-1/q2)
    cand <- rtrunc(1, spec="norm",a=2.1, b=100, mean = aw, sd = sqrt(bw))                                   
    alfa<-(exp( gHierar(cand,U))/exp( gHierar(last,U)))*(dtrunc(last,spec="norm",a=2.1, b=100, mean = aw, sd = sqrt(bw))/dtrunc(cand,spec="norm", a=2.1, b=100, mean = aw, sd = sqrt(bw)))
  }
  
  ifelse(runif(1) < min(alfa, 1), last <- cand, last<-last)
  return(last)
}

MHrhoCN <-function(last,U,cc,sigma2,nu,s0,s1,cens)
{
  last1 <- last/(1-last)
  desv <- .001
  cand <- rlnorm(1,meanlog=log(last1), sdlog=sqrt(desv))
  ver <- verosCN(U,cc,cens=cens,sigma2,nu)
  g <- function(r,ver)
  {
    r1 <- r/(1+r)
    ff<- (r1)^(s0-1)*(1-r1)^(s1-1)*ver*(1/(1+r)^2)
    return(ff)
  }
  alfa <- g(cand,ver)*cand/(g(last1,ver)*last1)
  ifelse(runif(1) < min(alfa, 1), last <- cand, last<-last)
  return(last)
}

BayCens <- function(start.iter,end.iter,j,cad, ...)
{
  pb <- txtProgressBar(start.iter,end.iter,
                       initial = start.iter, style=3, width=10,
                       char=ifelse((cad ==1||cad==3),"+","*"))
  Sys.sleep(0.5); setTxtProgressBar(pb, j)
  cat("\r")
  cat("\r")
}


