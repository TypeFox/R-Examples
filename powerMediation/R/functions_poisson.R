# created on May 7, 2015
# functions to calculate power/sample size for simple Poisson regression

# mu.T - mean exposure time 
# phi - a measure of overdispersion 
# 
# assume x1 ~ N(mu.x1, sigma2.x1)
sizePoisson = function(beta0, beta1, mu.x1, sigma2.x1, 
  mu.T=1, phi=1, 
  alpha=0.05, power=0.8)
{
  za=qnorm(1-alpha/2)
  zb=qnorm(power)
  denom=  mu.T*exp(beta0)*beta1^2

  part1 = beta1*mu.x1+beta1^2*sigma2.x1/2
  vb1.a= exp(-part1)/sigma2.x1
  vb1.0=1/sigma2.x1

  numer = phi*(za*sqrt(vb1.0)+zb*sqrt(vb1.a))^2

  N=ceiling(numer/denom)

  return(N)
}

powerPoisson = function(beta0, beta1, mu.x1, sigma2.x1, 
  mu.T=1, phi=1, 
  alpha=0.05, N=50)
{
  za=qnorm(1-alpha/2)
  part1 = beta1*mu.x1+beta1^2*sigma2.x1/2
  vb1.a= exp(-part1)/sigma2.x1
  vb1.0=1/sigma2.x1

  numer1=sqrt(N*mu.T*exp(beta0)*beta1^2/phi)
  numer2=za*sqrt(vb1.0)
  denom=sqrt(vb1.a)

  power=pnorm( (numer1-numer2)/denom )

  return(power)
}

## example
#exp.b0=1
## exp(b1)/exp(b0)
#eb1.eb0=1.3
#mu.T=1
#phi=1
#alpha=0.05
#mu.x1=3.2
#sigma2.x1=2.1^2
#
#print(sizePoisson(beta0=log(exp.b0), beta1=log(eb1.eb0*exp.b0), 
#  mu.x1=mu.x1, sigma2.x1=sigma2.x1, 
#  mu.T=mu.T, phi=phi, 
#  alpha=alpha, power=0.8))
#
#NVec=seq(from=5, to=50, by=5)
#len=length(NVec)
#cat("\nlen=", len, "\n")
#
#cat("\n******* exp(beta1)/exp(beta0) = 1.3 >>>>>>>>>>>>\n")
#powerVec=lapply(1:len, function(i){
#  power.i=powerPoisson(beta0=log(exp.b0), beta1=log(eb1.eb0*exp.b0), 
#    mu.x1=mu.x1, sigma2.x1=sigma2.x1, 
#    mu.T=mu.T, phi=phi, 
#    alpha=alpha, N=NVec[i])
#  return(power.i)
#})
#
#mat=cbind(powerVec, NVec)
#print(mat)
#
#eb1.eb0=1.5
#cat("\n******* exp(beta1)/exp(beta0) = 1.5 >>>>>>>>>>>>\n")
#powerVec=lapply(1:len, function(i){
#  power.i=powerPoisson(beta0=log(exp.b0), beta1=log(eb1.eb0*exp.b0), 
#    mu.x1=mu.x1, sigma2.x1=sigma2.x1, 
#    mu.T=mu.T, phi=phi, 
#    alpha=alpha, N=NVec[i])
#  return(power.i)
#})
#
#mat=cbind(powerVec, NVec)
#print(mat)
#
#
