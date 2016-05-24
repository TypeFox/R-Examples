
##########################################################################################
#PAUSE BETWEEN PLOTS
##########################################################################################

readkey <- function()
{
  cat ("Press [enter] to continue")
  line <- readline()
}

##########################################################################################
#SIMULATE DATA FROM THE PROPOSED MODEL (LMM with Normal random effect and ALD error term)
##########################################################################################

genY = function(nj,x,z,beta,sigmae,D,p)
{
  n = length(nj)
  N = sum(nj)
  d = dim(x)[2]
  q = dim(z)[2]
  z = matrix(z,N,q)  
  zb = c()
  
  for (j in 1:n)
  {
    z1 = matrix(z[(sum(nj[1:j-1])+1):(sum(nj[1:j])),],ncol=q)
    b  = rmvnorm(n=1,mean=rep(0,q),sigma=as.matrix(D))
    zb = rbind(zb,z1%*%t(b))
  }
  
  e  = rALD(n = N,mu = 0,sigma = sigmae,p = p)
  y  = x%*%beta + zb + e
  return(y)
}

##########################################################################################
#ASSYMETRIC LAPLACE DISTRIBUTION (d,p,q and r functions)
##########################################################################################

dALD = function(y,mu=0,sigma=1,p=0.5)
{
  if(length(y) == 0) stop("y must be provided.")
  if(sum(y[is.na(y)==TRUE]) > 0) stop("There are some NA's values in y.")
  if(sigma <= 0) stop("sigma must be a positive number.")
  if(p >= 1 | p <= 0) stop("p must be a real number in (0,1).")
  if(abs(mu) ==Inf) stop("mu must be a finite real number.")
  
  densi = ifelse(test=y<mu,yes=(p*(1-p)/sigma)*exp((1-p)*(y-mu)/sigma),no=(p*(1-p)/sigma)*exp(-(p)*(y-mu)/sigma))
  return(densi)
}

pALD = function(q,mu=0,sigma=1,p=0.5,lower.tail=TRUE)
{
  if(length(q) == 0) stop("q must be provided.")
  if(sigma <= 0) stop("sigma must be a positive number.")
  if(p >= 1 | p <= 0) stop("p must be a real number in (0,1).")
  if(abs(mu) ==Inf) stop("mu must be a finite real number.")
  if(lower.tail != TRUE && lower.tail != FALSE) stop("lower.tail must be TRUE or FALSE.")
  
  acum = ifelse(test=q<mu,yes=p*exp((1-p)*(q-mu)/sigma),no=1-(1-p)*exp(-(p)*(q-mu)/sigma))
  ifelse(test=lower.tail == TRUE,yes=return(acum),no=return(1-acum))
}

qALD = function(prob,mu=0,sigma=1,p=0.5,lower.tail=TRUE)
{
  if(length(prob) == 0) stop("prob must be provided.")
  if(sigma <= 0) stop("sigma must be a positive number.")
  if(p >= 1 | p <= 0) stop("p must be a real number in (0,1).")
  if(abs(mu) ==Inf) stop("mu must be a finite real number.")
  if(lower.tail != TRUE && lower.tail != FALSE) stop("lower.tail must be TRUE or FALSE.")
  if(sum(prob > 1 | prob < 0) > 0) stop("All elements of prob must be real numbers in [0,1].")
  
  q = sapply(X=prob,FUN=function(prob,mu,sigma,p){ifelse(test=lower.tail == TRUE,yes=prob,no=(1-prob));
                                                  ifelse(test=prob<p,yes=mu+((sigma*log(prob/p))/(1-p)),no=mu-((sigma*log((1-prob)/(1-p)))/p))},mu=mu,sigma=sigma,p=p)
  return(q)
}

rALD = function(n,mu=0,sigma=1,p=0.5)
{
  if(length(n) == 0) stop("The sample size n must be provided.")
  if(n <= 0 || n%%1 !=0) stop("The sample size n must be a positive integer value.")
  if(sigma <= 0) stop("sigma must be a positive number.")
  if(p >= 1 | p <= 0) stop("p must be a real number in (0,1).")
  if(abs(mu) ==Inf) stop("mu must be a finite real number.")
  
  u = runif(n)
  r = sapply(X=u,FUN=qALD,mu=mu,sigma=sigma,p=p)
  return(r)
}

##########################################################################################
#b_i|y_i DENSITY (proportional)
##########################################################################################

densbiv=function(j,bi,beta,sigmae,D,p,y1,x1,z1,q=q,nj=nj)
{
  bi = matrix(bi,nrow = q,ncol = 1)
  dprod = 1
  for(k in 1:nj[j])
  {
    dprod = dprod * dALD(y=as.numeric(y1[k]),mu=as.numeric(x1[k,]%*%beta + z1[k,]%*%bi),sigma=as.numeric(sigmae),p=p)
  }
  dens=as.numeric(dprod*dmvnorm(x=as.numeric(bi),mean=matrix(rep(0,q),q,1),sigma=D))
  return(dens)
}

##########################################################################################
#METROPOLIS HASTINGS ALGORITHM in order to generate from b_i|y_i
##########################################################################################

MHbi2 = function(j,M,y1,x1,z1,bi,bibi,d,q,p,nj,beta=beta,sigmae=sigmae,D=D)
{  
  E_bi = bi
  V_bi = bibi - E_bi%*%t(E_bi)
  V_bi = V_bi + 10^-10
  GEN = matrix(NA,nrow=q,ncol=(M+1))
  count = 1
  GEN[,1]=rmvnorm(n = 1,mean = E_bi,sigma=V_bi)
  
  while(count <= M)
  {
    cand = rmvnorm(n=1,mean=as.vector(GEN[,count]),sigma=V_bi)
    
    c1 = densbiv(j=j,bi=cand,beta=beta,sigmae=sigmae,D=D,p=p,y1,x1,z1,q,nj)*dmvnorm(x=GEN[,count],mean=as.vector(cand),sigma=V_bi)
    c2 = densbiv(j=j,bi=GEN[,count],beta=beta,sigmae=sigmae,D=D,p=p,y1,x1,z1,q,nj)*dmvnorm(x=as.vector(cand),mean=as.vector(GEN[,count]),sigma=V_bi)
    alfa = c1/c2
    
    if(is.nan(alfa)+0==1) {alfa=0.0001}
    if (runif(1) < min(alfa, 1))
    {
      count = count + 1
      GEN[,count] = cand
    }
  }
  return(GEN[,2:(M+1)])
}

##########################################################################################
#ELIMINATION MATRIX (Transforms vec(A) into vech(A))
##########################################################################################

MElim = function(q)
{
  Dtest = matrix(data = 0,ncol = q,nrow = q)
  Dtest[upper.tri(Dtest, diag = T)]=1
  vDtest = as.vector(Dtest)
  Elim = matrix(data = 0,nrow = (q*(1+q)/2),ncol = q^2)
  count=1
  for(i in 1:q^2)
  {
    if(vDtest[i]==1)
    {
      Elim[count,i]=1
      count = count +1
    }
  }
  return(Elim)
}

##########################################################################################
#ESTIMATED MARGINAL LOG-LIKELIHOOD (Using Importance Sampling)
##########################################################################################

logveroIS = function(beta,sigmae,D,y,x,z,nj,bi,bibi,MIS,n,d,q,p)
{
  logvero = 0
  bi = matrix(data = bi,nrow = n,ncol = q)
  
  for(j in 1:n)
  {
    y1=y[(sum(nj[1:j-1])+1):(sum(nj[1:j]))]#APPROVED
    x1=matrix(x[(sum(nj[1:j-1])+1):(sum(nj[1:j])),],ncol=d)#APPROVED
    z1=matrix(z[(sum(nj[1:j-1])+1):(sum(nj[1:j])),],ncol=q)#APPROVED
    
    if(q==1){bibij=bibi[j]}else{bibij=bibi[j,,]}
    
    E_bi = bi[j,]
    V_bi = bibij - E_bi%*%t(E_bi)
    Bgen = rmvnorm(n = MIS,mean = E_bi,sigma=V_bi)
    sum = 0
    
    for(l in 1:MIS)
    {
      multden = 1
      for(k in 1:nj[j])
      {
        multden = multden * (dALD(y=y1[k],mu=as.numeric(x1[k,]%*%beta + z1[k,]%*%Bgen[l,]),sigma=sigmae,p=p))
      }
      multt = multden*(dmvnorm(x = Bgen[l,],mean=rep(0,q),sigma=D)/dmvnorm(x = Bgen[l,],mean = E_bi,sigma=V_bi))
      sum = sum + multt
    }
    prom    = sum/MIS
    logvero = logvero + log(prom)
  }
  return(logvero)
}