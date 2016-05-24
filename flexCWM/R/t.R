#####################################
## Re-parameterized t distribution ##
#####################################

denst <- function(x,mu=0,sd=1,df=100)
{
  if (df>300) { return (dnorm(x,mu,sd))}
  else{
    num <- gamma((df+1)/2)
    den <- gamma(df/2)*sqrt(df*pi*sd)*(1+1/df*(x-mu)^2/sd)^((df+1)/2)
    return(num/den)
  }
}

probt <- function(q,mu=0,sigma=1,df=100) 
{
  #ploum <- function(y) denst(y, mu, sigma, df)
  res <- sapply(q, function(i) integrate(denst, lower=-Inf, upper=i, mu=mu, sd=sigma, df=df)$value)
  return(res)  
}

quantt <- function(p,mu=0,sigma=1,df=100) 
{
  
  n   <- length(p) 
  res <- numeric(n)
  for(i in 1:n){
    
    f <- function(par,p,mu=mu,sigma=sigma,df=df) 
      (probt(q=par,mu=mu,sigma=sigma,df=df)-p)^2
    
    #res[i] <- optim(par=qnorm(p=p[i],mean=mu,sd=sqrt(sigma)),fn=f,p=p[i],method="BFGS")$par 
    #res[i] <- optimize(f, lower = -Inf, upper = Inf, p=p[i], mu=mu, sigma=sigma, df=df)$minimum
    res[i] <- optim(par=qnorm(p=p[i],mean=mu,sd=sqrt(sigma)),fn=f,p=p[i],mu=mu,sigma=sigma,df=df,method="Nelder-Mead")$par 
    
  }
  res  
  
}

randt <- function(n,mu=0,sigma=1,df=100){
  
  values <- runif(n)
  sapply(1:n, function(i) quantt(p=values[i],mu=mu,sigma=sigma,df=df))
  
}
###################################
## Find df in the generic M-step ##
###################################

finddf <- function(z,v,dfold,mindf=2.001,maxdf=300)
{
  
  f <- function(par,z,v,dfold) 
    abs(-digamma(par/2)+log(par/2)+1+sum(z*(log(v)-v))/sum(z)+digamma((dfold+1)/2)-log((dfold+1)/2))
  
  df  <- optimize(f, c(mindf,maxdf), tol = 0.00001, z=z , v=v, dfold=dfold)
  
  return(df$minimum)  
}
