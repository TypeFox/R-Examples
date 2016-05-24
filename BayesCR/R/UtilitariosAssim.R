cdfSNI<-function(x,mu,sigma2,lambda,nu,type="SN")
{
  n <- length(x)
  resp<-matrix(0,n,1)
  if(type=="Normal")
  {
    resp <- pnorm((x-mu)/sqrt(sigma2))
    return(resp)
  }
  
  if(type=="T")
  {
    resp <- pt((x-mu)/sqrt(sigma2),df=nu)
    return(resp)
  }
  
  if(type=="SN")
  {
    delta <- lambda/sqrt(1+lambda^2)
    SIGMA <- matrix(c(sigma2,-delta*sqrt(sigma2),-delta*sqrt(sigma2),1),byrow=TRUE,ncol=2,nrow=2)
    if(length(mu)==1)
    {
      MU <- cbind(rep(mu,n),0)
    }
    if(length(mu)==n)
    {
      MU <- cbind(mu,0)
    }
    Y <- cbind(x,0)
    for(i in 1:n)
    {
      resp[i] <- 2*pmnorm(x=Y[i,],mean=MU[i,],varcov=SIGMA)
    }
    return(resp)
  }
  
  if(type=="ST")
  {
    delta <- lambda/sqrt(1+lambda^2)
    SIGMA <- matrix(c(sigma2,-delta*sqrt(sigma2),-delta*sqrt(sigma2),1),byrow=TRUE,ncol=2,nrow=2)
    if(length(mu)==1)
    {
      MU <- cbind(rep(mu,n),0)
    }
    if(length(mu)==n)
    {
      MU <- cbind(mu,0)
    }
    Y <- cbind(x,0)
    nu <- round(nu)
    for(i in 1:n)
    {
      resp[i] <- 2*pmt(x=Y[i,], mean = MU[i,], S=SIGMA, df=nu)
    }  
    return(resp)
    
  }
  if(type=="SSL")
  {
    cdf<- function(y)
    {
      f <- function(u) 2*nu*u^(nu - 1)*dnorm(y,mu,sqrt(u^(-1)*sigma2))*pnorm(u^(1/2)*lambda*(y-mu)/sqrt(sigma2))
      cdf <- integrate(Vectorize(f),0,1)$value
    }
    densidade <- as.numeric(cdf(x))
    resp<-as.numeric(integrate(Vectorize(cdf),-Inf,x)$value)
    return(list(pdf=densidade,cdf=resp))
  }
}

pdfSNI<-function(y,mu,sigma2,lambda,nu,type="SN")
{
  resp<-matrix(0,length(y),1)
  if(type=="Normal")
  {
    resp <- dnorm((y-mu)/sqrt(sigma2))/sqrt(sigma2)
  }
  
  if(type=="T")
  {
    resp <- dt((y-mu)/sqrt(sigma2),df=nu)/sqrt(sigma2)
  }
  
  if(type=="SN")
  {
    resp<-2*dnorm((y-mu)/sqrt(sigma2))*pnorm(lambda*(y-mu)/sqrt(sigma2))/sqrt(sigma2)
  }
  if(type=="ST")
  {
    z=(y-mu)/sqrt(sigma2)
    resp=2*dt(z,df=nu)*pt(sqrt(nu+1)*lambda*z/sqrt(nu+z^2),df=nu+1)/sqrt(sigma2)
  }
  if(type=="SSL")
  {
    f <- function(u){ 2*nu*u^(nu - 1)*dnorm(y,mu,sqrt(u^(-1)*sigma2))*pnorm(u^(1/2)*lambda*(y-mu)/sqrt(sigma2))}
    resp <- integrate(Vectorize(f),0,1)$value
  }

  return(resp)
}

dt.ls <- function(x, loc , sigma2 = 1,shape=1, nu = 4)
{
  d <- (x - loc)/sqrt(sigma2)
  dens <- 2*dt(d, df = nu)*pt(sqrt((1+nu)/(d^2+nu))*d*shape,1+nu)/sqrt(sigma2)
  return(dens)
}
                                       
MHnuST<-function(last,U,lambda)
{
	gHierar <- function(nu,U,lambda)
 	{
  	 	n <- length(U)
  		ff <- (-lambda*nu) + 0.5*n*nu*log(nu/2)+(0.5*nu)*sum(log(U)-U)-n*log(gamma(nu/2))
  	  return(ff)
 	}
      n <- length(U)
   	  Fonseca1<- deriv(~ (-lambda*nu) + 0.5*n*nu*log(nu/2)-n*log(gamma(nu/2)),c("nu"),function(nu){},hessian=TRUE)
   	  Fonseca2<- deriv(~(0.5*nu)*(log(U)-U),c("nu"),function(U,nu){},hessian=TRUE)
      aux1<- Fonseca1(last)
      aux2<- Fonseca2(U,last)                                                   
      q1<- attr(aux1,"gradient")[1]+sum(attr(aux2,"gradient"))
      q2<- attr(aux1,"hessian")[1]+sum(attr(aux2,"hessian"))
      aw<- last-q1/q2
      bw<- max(0.001,-1/q2)
      cand <- rtrunc(1, spec="norm",a=2.1, b=50, mean = aw, sd = sqrt(bw))                                   
      alfa<-(exp( gHierar(cand,U,lambda))/exp( gHierar(last,U,lambda)))*(dtrunc(last,spec="norm",a=2.1, b=50, mean = aw, sd = sqrt(bw))/dtrunc(cand,spec="norm", a=2.1, b=50, mean = aw, sd = sqrt(bw)))
      ifelse(runif(1) < min(alfa, 1), last <- cand, last<-last)
	    rejeitou <- 0
      ifelse(runif(1) < min(alfa, 1), rejeitou <- 0, rejeitou<-1)
	return(list(last=last,rejeitou=rejeitou))
}

rTN <- function(mu,sigma2,a,b){
	u <- runif(1)
	sd <- sqrt(sigma2)
	p <- u*(pnorm(b,mean=mu,sd=sd)-pnorm(a,mean=mu,sd=sd))+pnorm(a,mean=mu,sd=sd) 
	amostra<- qnorm(p,mean=mu,sd=sd)
	return(amostra)
}

rTGama <- function(alfa,beta,a,b){
	u <- runif(1)
	p <- u*(pgamma(b,shape=alfa,rate=beta)-pgamma(a,shape=alfa,rate=beta))+pgamma(a,shape=alfa,rate=beta) 
	amostra<- qgamma(p,shape=alfa,rate=beta)
	return(amostra)
}

rU <- function(alfa,beta){
	u <- runif(1)
	p <- u*pgamma(1,shape=alfa,rate=beta) 
	amostra<- qgamma(p,shape=alfa,rate=beta)
	return(amostra)
}

