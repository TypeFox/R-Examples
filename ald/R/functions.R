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


momentALD = function(k=1,mu=0,sigma=1,p=0.5)
{
  if(k <= 0 || k%%1 !=0) stop("The moment must be a positive integer value.")
  if(sigma <= 0) stop("sigma must be a positive number.")
  if(p >= 1 | p <= 0) stop("p must be a real number in (0,1).")
  if(abs(mu) ==Inf) stop("mu must be a finite real number.")
  
  
  moment = factorial(k)*(sigma^k)*p*(1-p)*(1/(p^(k+1)) + ((-1)^k)/((1-p)^(k+1)))
  return(moment)
}

meanALD = function(mu=0,sigma=1,p=0.5)
{
  mean = momentALD(k=1,mu=mu,sigma=sigma,p=p)
  return(mean)
}

varALD = function(mu=0,sigma=1,p=0.5)
{
  var = (momentALD(k=2,mu=mu,sigma=sigma,p=p) - (momentALD(k=1,mu=mu,sigma=sigma,p=p)^2))
  return(var)
}


skewALD = function(mu=0,sigma=1,p=0.5)
{
  skew = (momentALD(k=3,mu=mu,sigma=sigma,p=p) - 
            3*meanALD(mu,sigma,p)*varALD(mu,sigma,p) - 
            meanALD(mu,sigma,p)^3)/(varALD(mu,sigma,p)^(3/2))
  return(skew)
}

kurtALD = function(mu=0,sigma=1,p=0.5)
{
  kurt = (momentALD(k=4,mu=mu,sigma=sigma,p=p) - 
            4*meanALD(mu,sigma,p)*momentALD(k=3,mu=mu,sigma=sigma,p=p) + 
            6*(meanALD(mu,sigma,p)^2)*momentALD(k=2,mu=mu,sigma=sigma,p=p) -
            3*(meanALD(mu,sigma,p)^4))/(varALD(mu,sigma,p)^2)
  return(kurt)
}

absALD = function(sigma=1,p=0.5)
{
  if(sigma <= 0) stop("sigma must be a positive number.")
  if(p >= 1 | p <= 0) stop("p must be a real number in (0,1).")
  
  absmoment = (sigma)*((1 - 2*p + 2*p^2)/(p*(1-p)))
  return(absmoment)
}

# momentALD(k=3,mu=0,sigma=1,p=0.5)
# meanALD(mu=0,sigma=1,p=0.5)
# varALD(mu=0,sigma=1,p=0.5)
# skewALD(mu=0,sigma=1,p=0.5)
# kurtALD(mu=0,sigma=1,p=0.5)
# absALD(sigma=1,p=0.5)
# 
# likALD(y,mu=0,sigma=1,p=0.5,loglik = F)


likALD = function(y,mu=0,sigma=1,p=0.5,loglik=TRUE)
{
  if(loglik != TRUE && loglik != FALSE) stop("log must be TRUE or FALSE.")
  ifelse(test=loglik == FALSE,yes=return(prod(dALD(y,mu,sigma,p=p))),no=return(sum(log(dALD(y,mu,sigma,p=p)))))
}


minfunc = function(y,mu,p)
{
  if(length(y) == 0) stop("y must be provided.")
  if(sum(y[is.na(y)==TRUE]) > 0) stop("There are some NA's values in y.")
  if(p >= 1 | p <= 0) stop("p must be a real number in (0,1).")
  if(abs(mu) ==Inf) stop("mu must be a finite real number.")
  
  dif = y - mu
  return((sum(p*dif[dif>0]) - sum((1-p)*dif[dif<0])))
}

mleALD = function(y,initial=NA)
{
  if(is.na(initial) == TRUE)
  {
    p = 0.5
    OPT = optimize(f = minfunc,interval = c(-10000,10000),y=y,p=p)
    mu = OPT$minimum
    sigma = (1/length(y))*OPT$objective
    
  }else
  {
    if(abs(initial[1]) ==Inf) stop("mu must be a finite real number.")
    mu = initial[1]
    if(initial[2] <= 0) stop("sigma must be a positive number.")
    sigma = initial[2]
    if(initial[3] >= 1 | initial[3] <= 0) stop("p must be a real number in (0,1).")
    p = initial[3]
  }

  teta0 = c(mu,sigma,p)
  criterio = 1
  iter = 0
  while(max(criterio) > 10^-5)
  {
    p = optimize(f = likALD,interval = c(0.001,0.999),y=y,mu=mu,sigma=sigma,loglik=T,maximum = T)$maximum
    OPT = optimize(f = minfunc,interval = c(-10000,10000),y=y,p=p)
    mu = OPT$minimum
    sigma = (1/length(y))*OPT$objective
    tetaT = c(mu,sigma,p)
    criterio = (teta0-tetaT)^2
    teta0 = tetaT
    iter = iter +1
  }
  
  return(list(iter = iter, par = teta0))
}

# y = rALD(10000,mu = -323,sigma = 40,p = 0.9)
# res = mleALD(y)
# seqq = seq(-3000,0)
# dens = dALD(y=seqq,mu=res$par[1],sigma=res$par[2],p=res$par[3])
# hist(y,breaks=50,freq = F,ylim=c(0,max(dens)))
# lines(seqq,dens,type="l",lwd=2,col="red",xlab="x",ylab="f(x)", main="ALD Density function")

