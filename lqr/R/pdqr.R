##########################################################################################
# (d,p,q and r functions)
##########################################################################################

dSKD = function(y,mu=0,sigma=1,p=0.5,dist = "normal",nu="",gama="")
{
  if(any(y==Inf)) stop("y must be a finite real number.")
  if(sigma <= 0) stop("sigma must be a positive number.")
  if(p >= 1 | p <= 0) stop("p must be a real number in (0,1).")
  if(abs(mu) ==Inf) stop("mu must be a finite real number.")
  if(gama != "" && (gama >= 1 | gama <= 0) && dist == "cont") stop("nu must be a real number in (0,1)")
  if(nu != "" && (nu >= 1 | nu <= 0) && dist == "cont") stop("nu must be a real number in (0,1)")
  #if(nu != "" && (nu >= 100 | nu < 2) && dist == "t") stop("nu must be a positive real number at least 2.")
  if(nu != "" && (nu < 2) && dist == "t") stop("nu must be a positive real number at least 2.")
  if(nu != "" && nu <= 0 && dist == "slash") stop("nu must be a positive real number.")
  if(dist != "" && dist != "normal" && dist != "t" && dist != "laplace" && dist != "slash" && dist != "cont") stop("The dist values are normal, t, laplace, slash or cont.")
  
  #Default values
  if(dist == ""){dist = "normal"}
  if(nu=="" && dist == "t"){nu=4}
  if(nu=="" && dist == "slash"){nu=2}
  if(nu=="" && dist == "cont"){nu=0.1}
  if(gama=="" && dist == "cont"){gama=0.1}
  
  if(dist == "normal"){return(densN(y,mu = mu,sigma = sigma,p = p))}
  if(dist == "t"){return(densT(x = y,mu = mu,sigma = sigma,nu = nu,p = p))}
  if(dist == "laplace"){return(densL(y,mu = mu,sigma = sigma,p = p))}
  if(dist == "slash"){return(densSl(y,mu = mu,sigma = sigma,nu = nu,p = p))}
  if(dist == "cont"){return(densNC(y,mu = mu,sigma = sigma,nu = nu,gama = gama,p = p))}
}

# seqq = seq(from=-20,to = 10,length.out = 1000)
# vals = dSKD(x=seqq,mu=0,sigma=1,p=0.85,dist="laplace")
# plot(seqq,vals,type="l")

pSKD = function(q,mu=0,sigma=1,p=0.5,dist = "normal",nu="",gama="",lower.tail=TRUE)
{
  #if(any(q==Inf)) stop("q must be a finite real number.")
  if(sigma <= 0) stop("sigma must be a positive number.")
  if(p >= 1 | p <= 0) stop("p must be a real number in (0,1).")
  if(abs(mu) ==Inf) stop("mu must be a finite real number.")
  if(gama != "" && (gama >= 1 | gama <= 0) && dist == "cont") stop("nu must be a real number in (0,1)")
  if(nu != "" && (nu >= 1 | nu <= 0) && dist == "cont") stop("nu must be a real number in (0,1)")
  if(nu != "" && (nu >= 100 | nu < 2) && dist == "t") stop("nu must be a positive real number at least 2.")
  if(nu != "" && nu <= 0 && dist == "slash") stop("nu must be a positive real number.")
  if(dist != "" && dist != "normal" && dist != "t" && dist != "laplace" && dist != "slash" && dist != "cont") stop("The dist values are normal, t, laplace, slash or cont.")
  if(lower.tail != TRUE && lower.tail != FALSE) stop("lower.tail must be TRUE or FALSE.")
  
  #Default values
  if(dist == ""){dist = "normal"}
  if(nu=="" && dist == "t"){nu=4}
  if(nu=="" && dist == "slash"){nu=2}
  if(nu=="" && dist == "cont"){nu=0.1}
  if(gama=="" && dist == "cont"){gama=0.1}
  
  #if(dist == "normal"){out = integrate(f = densN,lower = -Inf,upper = x,mu = mu,sigma = sigma,p = p)$value}
  if(dist == "normal"){out = as.double(as.matrix(sapply(X = q,FUN = integrate,f = densN,lower = -Inf,mu = mu,sigma = sigma,p = p))[1,])}
  if(dist == "t"){out = as.double(as.matrix(sapply(X = q,FUN = integrate,f = densT,lower = -Inf,mu = mu,sigma = sigma,p = p,nu=nu))[1,])}
  if(dist == "laplace"){out = as.double(as.matrix(sapply(X = q,FUN = integrate,f = densL,lower = -Inf,mu = mu,sigma = sigma,p = p))[1,])}
  if(dist == "slash"){out = as.double(as.matrix(sapply(X = q,FUN = integrate,f = densSl,lower = -Inf,mu = mu,sigma = sigma,p = p,nu=nu))[1,])}
  if(dist == "cont"){out = as.double(as.matrix(sapply(X = q,FUN = integrate,f = densNC,lower = -Inf,mu = mu,sigma = sigma,p = p,nu=nu,gama=gama))[1,])}
  ifelse(test=lower.tail == TRUE,yes=return(out),no=return(1-out))
}

# seqq = seq(from=-20,to = 10,length.out = 1000)
# vals = pSKD(x=seqq,mu=0,sigma=1,p=0.75,dist="t")
# plot(seqq,vals,type="l")


qSKD = function(prob,mu=0,sigma=1,p=0.5,dist = "normal",nu="",gama="",lower.tail=TRUE)
{
  if(any(prob > 1) | any(prob < 0)) stop("prob must be a real number in (0,1).")
  if(sigma <= 0) stop("sigma must be a positive number.")
  if(p > 1 | p < 0) stop("p must be a real number in [0,1].")
  if(abs(mu) ==Inf) stop("mu must be a finite real number.")
  if(gama != "" && (gama >= 1 | gama <= 0) && dist == "cont") stop("nu must be a real number in (0,1)")
  if(nu != "" && (nu >= 1 | nu <= 0) && dist == "cont") stop("nu must be a real number in (0,1)")
  if(nu != "" && (nu >= 100 | nu < 2) && dist == "t") stop("nu must be a positive real number at least 2.")
  if(nu != "" && nu <= 0 && dist == "slash") stop("nu must be a positive real number.")
  if(dist != "" && dist != "normal" && dist != "t" && dist != "laplace" && dist != "slash" && dist != "cont") stop("The dist values are normal, t, laplace, slash or cont.")
  if(lower.tail != TRUE && lower.tail != FALSE) stop("lower.tail must be TRUE or FALSE.")
  
  #Default values
  if(dist == ""){dist = "normal"}
  if(nu=="" && dist == "t"){nu=4}
  if(nu=="" && dist == "slash"){nu=2}
  if(nu=="" && dist == "cont"){nu=0.1}
  if(gama=="" && dist == "cont"){gama=0.1}
  if(lower.tail == FALSE){prob = 1 - prob}
  
  mu0 = mu
  mu  = 0
  
  f <- function(x) pSKD(x,mu=mu,sigma=sigma,p=p,dist = dist,nu = nu,gama = gama)
  f.inv <- inverse(f,lower=-Inf,upper=Inf)
  
  kind  = rep(NA,length(prob)) 
  tempk = rep(NA,length(prob))
  
  for(k in 1:length(prob))
  {
    temp = tryCatch(f.inv(prob[k]),error=function(e) e)
    if((inherits(temp,"error")))
    {
      kind[k] = k;
    }
    else{tempk[k] = temp}
  }
  
  #print(length(kind[!is.na(kind)]))
  
  tempk[kind[!is.na(kind)]] = quantile(x = rSKD(n=10^6,mu=mu,sigma=sigma,p=p,dist=dist,nu=nu),probs = prob[kind[!is.na(kind)]])
  
  return(round(as.double(tempk),digits = 5)+mu0)
}

# seqq = seq(from=0.001,to = 0.999,length.out = 100)
# idf = qSKD(prob=seqq,mu=50,sigma=1,p=0.75,dist="normal")
# plot(seqq,idf,type="l")


rSKD = function(n,mu=0,sigma=1,p=0.5,dist = "normal",nu="",gama="")
{
  if(length(n) == 0) stop("The sample size n must be provided.")
  if(n <= 0 || n%%1 !=0) stop("The sample size n must be a positive integer value.")
  if(sigma <= 0) stop("sigma must be a positive number.")
  if(p >= 1 | p <= 0) stop("p must be a real number in (0,1).")
  if(abs(mu) ==Inf) stop("mu must be a finite real number.")
  if(gama != "" && (gama >= 1 | gama <= 0) && dist == "cont") stop("nu must be a real number in (0,1)")
  if(nu != "" && (nu >= 1 | nu <= 0) && dist == "cont") stop("nu must be a real number in (0,1)")
  if(nu != "" && (nu >= 100 | nu < 2) && dist == "t") stop("nu must be a positive real number at least 2.")
  if(nu != "" && nu <= 0 && dist == "slash") stop("nu must be a positive real number.")
  if(dist != "" && dist != "normal" && dist != "t" && dist != "laplace" && dist != "slash" && dist != "cont") stop("The dist values are normal, t, laplace, slash or cont.")
  
  return(genSLK(n,mu,sigma,p,dist,nu,gama))
}

# x = rSKD(n=100000,mu=0,sigma=2,p=0.2,dist="cont",gama=0.3)
# hist(x,breaks = 100)


#package.skeleton(list = c("dSKD","pSKD","qSKD","rSKD"), name = "mypkg2")
