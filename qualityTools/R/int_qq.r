

.confintbeta= function(thethas, varmatrix, alpha) {

    th1= thethas[[1]]
    th2 =  thethas[[2]]

    prozent=c(0.001, seq(0.01,0.09, 0.01 ), seq(0.1,0.9,0.1), seq(0.91, 0.99, 0.01), 0.999)

    perzentile=qbeta(prozent, th1, th2)

    h=1e-6
    dFdth1=(qbeta(prozent, th1, th2)-qbeta(prozent, th1+h, th2))/h
    dFdth2=(qbeta(prozent, th1, th2)-qbeta(prozent, th1, th2+h))/h

    Var = varmatrix[[1, 1]]*dFdth1^2 + 2*varmatrix[[1, 2]]*dFdth1*dFdth2 + varmatrix[[2, 2]]*dFdth2^2
    zalpha=qnorm(1-alpha/2)
    halfwidth = zalpha*sqrt(Var)

    lci=perzentile-halfwidth
    uci=perzentile+halfwidth

	bounds = list(lci, uci, perzentile)


return (bounds)
}

.confintcauchy = function(thethas, varmatrix, alpha) {

    th1= thethas[[1]]
    th2 =  thethas[[2]]

    prozent=c(0.001, seq(0.01,0.09, 0.01 ), seq(0.1,0.9,0.1), seq(0.91, 0.99, 0.01), 0.999)

    perzentile=qcauchy(prozent, th1, th2)

    h=1e-6
    dFdth1=(qcauchy(prozent, th1, th2)-qcauchy(prozent, th1+h, th2))/h
    dFdth2=(qcauchy(prozent, th1, th2)-qcauchy(prozent, th1, th2+h))/h

    Var = varmatrix[[1, 1]]*dFdth1^2 + 2*varmatrix[[1, 2]]*dFdth1*dFdth2 + varmatrix[[2, 2]]*dFdth2^2
    zalpha=qnorm(1-alpha/2)
    halfwidth = zalpha*sqrt(Var)

    lci=perzentile-halfwidth
    uci=perzentile+halfwidth

	bounds = list(lci, uci, perzentile)


return (bounds)
}

.confintexp=function(thethas, varmatrix, alpha) {
  lambda=thethas[[1]]
  prozent=c(0.001, seq(0.01,0.09, 0.01 ), seq(0.1,0.9,0.1), seq(0.91, 0.99, 0.01), 0.999)
  perzentile=qexp(prozent, lambda)
  logPerzentile = log(perzentile)
   zalpha=qnorm(1-alpha/2)
  halfwidth = zalpha*sqrt(varmatrix[[1, 1]]/lambda^2)
   lci = exp(logPerzentile - halfwidth);
   uci = exp(logPerzentile + halfwidth);

   bounds = list(lci, uci, perzentile)

return (bounds)

}

.confintgamma= function(thethas, varmatrix, alpha) {
  th1= thethas[[1]]
  th2 =  thethas[[2]]
  prozent=c(0.001, seq(0.01,0.09, 0.01 ), seq(0.1,0.9,0.1), seq(0.91, 0.99, 0.01), 0.999)
  perzentile=qgamma(prozent, th1, th2)

    h=1e-6
    dFdth1=(qgamma(prozent, th1, th2)-qgamma(prozent, th1+h, th2))/h
    dFdth2=(qgamma(prozent, th1, th2)-qgamma(prozent, th1, th2+h))/h

    Var = varmatrix[[1, 1]]*dFdth1^2 + 2*varmatrix[[1, 2]]*dFdth1*dFdth2 + varmatrix[[2, 2]]*dFdth2^2
    zalpha=qnorm(1-alpha/2)
    halfwidth = zalpha*sqrt(Var)

    lci=perzentile-halfwidth
    uci=perzentile+halfwidth

  bounds = list(lci, uci, perzentile)

return (bounds)
}

.confintlnorm=function(thethas, varmatrix, alpha){
    th1= thethas[[1]]
    th2 =  thethas[[2]]

prozent=c(0.001, seq(0.01,0.09, 0.01 ), seq(0.1,0.9,0.1), seq(0.91, 0.99, 0.01), 0.999)
perzentile=qlnorm(prozent, th1, th2)

zp=qnorm(prozent)

varPerzentile = varmatrix[[1, 1]]+2*varmatrix[[1, 2]]*zp+varmatrix[[2, 2]]*zp*zp

   zalpha=qnorm(1-alpha/2)
lci=log(perzentile)-zalpha*sqrt(varPerzentile)
uci=log(perzentile)+zalpha*sqrt(varPerzentile)

bounds = list(exp(lci), exp(uci), perzentile)

return (bounds)

}

.confintlogis= function(thethas, varmatrix, alpha) {

    th1= thethas[[1]]
    th2 =  thethas[[2]]

    prozent=c(0.001, seq(0.01,0.09, 0.01 ), seq(0.1,0.9,0.1), seq(0.91, 0.99, 0.01), 0.999)

    perzentile=qlogis(prozent, th1, th2)

    h=1e-6
    dFdth1=(qlogis(prozent, th1, th2)-qlogis(prozent, th1+h, th2))/h
    dFdth2=(qlogis(prozent, th1, th2)-qlogis(prozent, th1, th2+h))/h

    Var = varmatrix[[1, 1]]*dFdth1^2 + 2*varmatrix[[1, 2]]*dFdth1*dFdth2 + varmatrix[[2, 2]]*dFdth2^2
    zalpha=qnorm(1-alpha/2)
    halfwidth = zalpha*sqrt(Var)

    lci=perzentile-halfwidth
    uci=perzentile+halfwidth

	bounds = list(lci, uci, perzentile)


return (bounds)
}

.confintnorm=function(thethas, varmatrix, alpha){

    prozent=c(0.001, seq(0.01,0.09, 0.01 ), seq(0.1,0.9,0.1), seq(0.91, 0.99, 0.01), 0.999)
zp=qnorm(prozent)
perzentile=qnorm(prozent, thethas[[1]], thethas[[2]])

varPerzentile = varmatrix[[1, 1]]+2*varmatrix[[1, 2]]*zp+varmatrix[[2, 2]]*zp*zp

zalpha=qnorm(1-alpha/2)
lci=perzentile-zalpha*sqrt(varPerzentile)
uci=perzentile+zalpha*sqrt(varPerzentile)

bounds = list(lci, uci, perzentile)


return (bounds)
}

.confintweibull= function(thethas, varmatrix, alpha) {
    th1= thethas[[1]]
    th2 =  thethas[[2]]

    prozent=c(0.001, seq(0.01,0.09, 0.01 ), seq(0.1,0.9,0.1), seq(0.91, 0.99, 0.01), 0.999)
perzentile=qweibull(prozent, th1, th2)
q=-log(1-prozent)
logPerzentile=log(perzentile)
logq=log(q)
dB=1/th2
dA=-1/(th1^2)

Var = varmatrix[[1, 1]]*(dA*logq)^2 + 2*varmatrix[[1, 2]]*dB*dA*logq + varmatrix[[2, 2]]*dB^2
zalpha=qnorm(1-alpha/2)
halfwidth = zalpha*sqrt(Var)


lci=exp(logPerzentile-halfwidth)
uci=exp(logPerzentile+halfwidth)

bounds = list(lci, uci, perzentile)

# print(data.frame(prozent, uci, perzentile, lci))

return (bounds)
}

.gamma3 = function(data) {
     n=length(data)
     data=sort(data)

     pEmp= (seq(1:n)-0.5)/n

     weight = 1 / sqrt(pEmp*(1-pEmp))

     thld = .99*min(data)
     shape=1
     scale=1

   gammaEst = function(param) {
    return( sum(weight*(pgamma(data-param[3], shape = exp(param[1]), scale = exp(param[2]))-pEmp)^2) )
  }

     paramEst = optim(c(shape, scale, thld), gammaEst, method = "Nelder-Mead")
     paramEst = paramEst$par
     return(list(shape = exp(paramEst[1]), scale = exp(paramEst[2]), threshold = paramEst[3]))
}

.lognormal3 = function(data) {

  n=length(data)
  data=sort(data)
  #compute the empirical cumulative distribution function of the data
  pEmp= (seq(1:n)-0.5)/n
   # will minimize the weighted sum of squared distances
   # so compute weights
  weight = 1 / sqrt(pEmp*(1-pEmp))

  # initial values for optimization
  thld = .99*min(data)
  mu0 = mean(log(data-thld))
  sigma0 = sd(log(data-thld))


  lnEst = function(param) {
    return( sum(weight*(plnorm(data-param[3], meanlog = param[1], sdlog = exp(param[2]))-pEmp)^2) )
  }

  logSigma0=log(sigma0)
  # optimize gammaEst using optim function
  paramEst = optim(c(mu0,logSigma0, thld), lnEst, method = "Nelder-Mead")
  param = paramEst$par

  return(list(meanlog = param[1], sdlog = exp(param[2]), threshold = param[3]))
}

.weibull3 = function(x)
{
#  if(any(x < 0))                                                               ####
#    stop("x must be positive")                                                 ####

  n = length(x)
  x = sort(x)
  p = ((1:n)-0.5)/n
  interval = c(0.75*min(x), 0.9999*min(x))

  wb3RSquared = function(th)
  {
      return(summary(lm(log(x-th) ~ log(-log(1-p))))$r.squared)
  }

  th = (optimize(wb3RSquared, interval = interval, maximum = TRUE))$maximum

  lm.1 = lm(log(x-th) ~ log(-log(1-p)))
  estimates = list(shape = 1/coef(lm.1)[[2]], scale = exp(coef(lm.1)[[1]]), threshold = th)
  return(estimates)
}