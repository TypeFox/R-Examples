# Models (Sobel, MacKinnon et al.):
# mi=theta0+theta1*xi+ei, ei~N(0, sigma.e^2)
# yi=gamma+lambda*mi+lambda2*xi+epsiloni, epsiloni~N(0, sigma.epsilone^2)
# H0: theta_1=0 and lambda=0 versus theta_1 neq 0 or lambda neq 0
#
# Model (Vittinghoff and McCulloch, 2009)
# (1) y=a0+a1*x+epsilon
# (2) y=b0+b1*x+b2*m+e, e~N(0, sigma.e^2)
# H0: b2=0 vs Ha: b2 neq 0
# 
#
# p-value and CI for mediation effect based on Sobel's test
testMediation.Sobel <- function(theta.1.hat, lambda.hat, 
  sigma.theta1, sigma.lambda, alpha=0.05)
{
  part1<-(theta.1.hat*sigma.lambda)^2
  part2<-(lambda.hat*sigma.theta1)^2

  sigma.hat<-sqrt(part1+part2)
  Z.obs<-theta.1.hat*lambda.hat/sigma.hat

  pval<-2*(1-pnorm(abs(Z.obs)))

  za2<-qnorm(1-alpha/2)

  t1<-theta.1.hat*lambda.hat
  t2<-za2*sigma.hat

  CI.low<-t1-t2
  CI.upp<-t1+t2

  res<-list(pval=pval, CI.low=CI.low, CI.upp=CI.upp)
  return(res)
}

# p-value and CI for mediation effect based on MacKinnon's test
testMediation.MacKinnon <- function(n, theta.1.hat, lambda.hat, 
  sigma.theta1.hat, sigma.lambda.hat, alpha=0.05)
{
  numer<-theta.1.hat*lambda.hat
  denom<-sigma.theta1.hat*sigma.lambda.hat
  C.obs<-numer/denom

  Z1<-theta.1.hat/sigma.theta1.hat
  Z2<-lambda.hat/sigma.lambda.hat
  tmp<-CDFprod(z=abs(C.obs), xmux=0, xmuy=0, rho=0, verbose=FALSE)$prob
  pval<-2*(1-tmp)

  c1<-qprod(p=alpha/2, xmux=Z1, xmuy=Z2, rho=0, verbose=FALSE)$z-C.obs
  c2<-qprod(p=1-alpha/2, xmux=Z1, xmuy=Z2, rho=0, verbose=FALSE)$z-C.obs

  CI.low<-numer-c2*denom
  CI.upp<-numer-c1*denom

  res<-list(pval=pval, CI.low=CI.low, CI.upp=CI.upp)
  return(res)
}

######################################

# based on MacKinnon et al.'s Test 
# mi = theta0 + theta1*xi + ei, ei~N(0, sigma.e^2)
# yi = gamma + lambda * mi + lambda2*xi+epsiloni, epsiloni~N(0, sigma.epsilon^2)
# rho.mx=cor(x, m), var(m)=sigma.m^2, var(x)=sigma.x^2
#
# H0: theta1*lambda = 0 vs Ha: theta1*lambda=theta.1a*lambda.a neq 0
powerMediation.MacKinnon <- function(n, theta.1a, lambda.a, sigma.x, 
  sigma.m, rho2.mx, sigma.e, sigma.epsilon, alpha=0.05, verbose=TRUE)
{
    if(rho2.mx >= 1 || rho2.mx < 0)
    {
        stop("rho2.mx should be in the range [0, 1)")
    }

    numer <- n * theta.1a * lambda.a * sigma.x * sigma.m * sqrt(1-rho2.mx)
    denom <- sigma.e * sigma.epsilon
    delta <- numer / denom

    # expection of Z1 and Z2
    EZ1=theta.1a*sigma.x*sqrt(n)/sigma.e
    EZ2=lambda.a*sigma.m*sqrt(n*(1-rho2.mx))/sigma.epsilon
  
    alpha2<-alpha/2
    ca2 <- qprod(p = 1 - alpha2, xmux = 0, xmuy = 0, rho = 0, verbose = FALSE)
    ca2 <- ca2$z

    p1 <- CDFprod(z=ca2-delta, xmux = EZ1, xmuy = EZ2, rho = 0, verbose = FALSE) 
    p1 <- p1$prob
    p2 <- CDFprod(z=-ca2-delta, xmux = EZ1, xmuy = EZ2, rho = 0, 
      verbose = FALSE) 
    p2 <- p2$prob
    
    power<-1-p1+p2
  
    res<-list(power=power, delta=delta)
    if(verbose)
    { 
        print(res$power) 
    }
    invisible(res)
}



# based on Sobel First-Order Test 
# mi = theta0 + theta1*xi + ei, ei~N(0, sigma.e^2)
# yi = gamma + lambda * mi + epsiloni, epsiloni~N(0, sigma.epsilon^2)
# H0: theta1*lambda = 0 vs Ha: theta1*lambda=theta.1a*lambda.a neq 0
powerMediation.Sobel <- function(n, theta.1a, lambda.a, sigma.x, 
  sigma.m, rho2.mx, sigma.epsilon, alpha=0.05, verbose=TRUE)
{
    if(rho2.mx >= 1 || rho2.mx < 0)
    {
        stop("rho2.mx should be in the range [0, 1)")
    }

    numer<-sqrt(n)*theta.1a*lambda.a
    # sigma.e^2=sigma.m^2*(1-rho.mx^2)
    denom<-sqrt(theta.1a^2*sigma.epsilon^2/(sigma.m^2*(1-rho2.mx))+lambda.a^2*sigma.m^2*(1-rho2.mx)/sigma.x^2)
    delta<-numer/denom
  
    alpha2<-alpha/2
    za2<-qnorm(1-alpha2)
    power<-1-pnorm(za2-delta)+pnorm(-za2-delta)
  
    res<-list(power=power, delta=delta)
    if(verbose)
    { 
        print(res$power) 
    }
    invisible(res)
}


# cumulative distribution of the product of two normal random variables 
# that bivariate normal distribution
# wrapper for Alan Miller's fortran code fnprod
# http://jblevins.org/mirror/amiller/
# pr(X*Y < z)
CDFprod <- function(z, xmux, xmuy, rho, verbose = TRUE)
{
    if(rho > 1 || rho < -1)
    {
        stop("rho should be in the range [-1, 1]")
    }
  
    res <- .Fortran("fnprod", 
         as.double(xmux), 
         as.double(xmuy), 
         as.double(rho),
         as.double(z),
         answer = as.double(0),
         ier = as.integer(0), 
         abserr = as.double(0),
         last = as.integer(0),
         PACKAGE = "powerMediation")  
       
     # prob = Pr(X*Y < z), where (X, Y) has bivariate normal distr.
     resList <- list(prob = res$answer, 
         ier = res$ier,
         abserr = res$abserr, 
         last = res$last) 

    if(verbose)
    {
        print(resList$prob)
    }
    return(resList)
}

tmpCDFprod <- function(z, p, xmux, xmuy, rho)
{
    tmp <- CDFprod(z, xmux, xmuy, rho, verbose = FALSE)
    res <- tmp$prob - p
    return(res)
}
 
# quantile of the product of two normal random variables 
# that bivariate normal distribution.
# find z such that pr(X*Y < z)=p
qprod <- function(p,  xmux, xmuy, rho, z.lower=-1.0e+30, z.upper=1.0e+30,
  verbose = TRUE)
{
    if(p > 1 || p < 0)
    {
        stop("p should be in the range [0, 1]")
    }
    if(rho > 1 || rho < -1)
    {
        stop("rho should be in the range [-1, 1]")
    }
    if(z.lower >= z.upper)
    {
      stop("z.lower must be < z.upper")
    }

    res.uniroot <- uniroot(f = tmpCDFprod,
        interval=c(z.lower, z.upper),
        p = p, xmux = xmux, xmuy = xmuy, rho = rho)
  
    z = res.uniroot$root
    if(verbose)
    {
        print(z)
    }
  
    res<-list(z = z, res.uniroot = res.uniroot)
  
    invisible(res)
}


tmpSS.mediation.Sobel<-function(n, power, theta.1a, lambda.a, 
   sigma.x, sigma.m, rho2.mx, sigma.epsilon,
   alpha=0.05, verbose=FALSE)
{
  tmppower<-powerMediation.Sobel(n=n, theta.1a=theta.1a, 
    lambda.a=lambda.a, sigma.x=sigma.x, sigma.m=sigma.m,
    rho2.mx=rho2.mx, sigma.epsilon=sigma.epsilon,
    alpha=alpha, verbose=verbose)
  res<- tmppower$power-power
  return(res)
}

#####################################################
ssMediation.Sobel <- function(power, theta.1a, lambda.a, 
  sigma.x, sigma.m, rho2.mx, sigma.epsilon, 
  n.lower=1, n.upper=1.0e+30, 
  alpha = 0.05, verbose=TRUE)
{
  if(n.lower< 1)
  {
    stop("n.lower must be >= 1")
  }
  if(n.lower >= n.upper)
  {
    stop("n.lower must be < n.upper")
  }
  if(power<0 || power > 1)
  {
    stop("power must be in the range [0, 1]")
  }  
  if(rho2.mx >= 1 || rho2.mx < 0)
  {
      stop("rho2.mx should be in the range [0, 1)")
  }


  res.uniroot<-uniroot(f=tmpSS.mediation.Sobel, 
      interval=c(n.lower, n.upper),
      power=power, theta.1a=theta.1a, lambda.a=lambda.a, 
      sigma.x=sigma.x, sigma.m=sigma.m, rho2.mx=rho2.mx,
      sigma.epsilon=sigma.epsilon,
      alpha=alpha, verbose=FALSE)

  n.numeric<-res.uniroot$root

  res<-list(n=n.numeric, res.uniroot=res.uniroot)
  if(verbose)
  { print(res$n) }
  invisible(res)
}

tmpSS.mediation.MacKinnon<-function(n, power, theta.1a, lambda.a, 
   sigma.x, sigma.m, rho2.mx, sigma.e, sigma.epsilon,
   alpha=0.05, verbose=FALSE)
{
  tmppower<-powerMediation.MacKinnon(n=n, theta.1a=theta.1a, 
    lambda.a=lambda.a, sigma.x=sigma.x, sigma.m=sigma.m,
    rho2.mx=rho2.mx, sigma.e=sigma.e, sigma.epsilon=sigma.epsilon,
    alpha=alpha, verbose=verbose)
  res<- tmppower$power-power
  return(res)
}

#####################################################
ssMediation.MacKinnon <- function(power, theta.1a, lambda.a, 
  sigma.x, sigma.m, rho2.mx, sigma.e, sigma.epsilon, n.lower=1, n.upper=1e+30, 
  alpha = 0.05, verbose=TRUE)
{
  if(n.lower< 1)
  {
    stop("n.lower must be >= 1")
  }
  if(n.lower >= n.upper)
  {
    stop("n.lower must be < n.upper")
  }
  if(power<0 || power > 1)
  {
    stop("power must be in the range [0, 1]")
  }  
  if(rho2.mx >= 1 || rho2.mx < 0)
  {
      stop("rho2.mx should be in the range [0, 1)")
  }

  res.uniroot<-uniroot(f=tmpSS.mediation.MacKinnon, 
      interval=c(n.lower, n.upper),
      power=power, theta.1a=theta.1a, lambda.a=lambda.a, 
      sigma.x=sigma.x, sigma.m=sigma.m, rho2.mx=rho2.mx,
      sigma.e=sigma.e, sigma.epsilon=sigma.epsilon, 
      alpha=alpha, verbose=FALSE)

  n.numeric<-res.uniroot$root

  res<-list(n=n.numeric, res.uniroot=res.uniroot)
  if(verbose)
  { print(res$n) }
  invisible(res)
}

#########################################################
# for simple linear regression y=theta + lambda * x + epsilon
#########################################################

power.SLR<-function(n, lambda.a, sigma.x, sigma.y, alpha=0.05, verbose=TRUE)
{
  if(n <= 2)
  {
    stop("n must be > 2")
  }
  if(lambda.a > sigma.y/sigma.x || lambda.a < -sigma.y/sigma.x)
  {
    stop("lambda.a must be in the range [-sigma.y/sigma.x, sigma.y/sigma.x]")
  }

  s=sqrt(sigma.y^2-(lambda.a*sigma.x)^2)
  delta = lambda.a*sigma.x*sqrt(n)/s

  rho=lambda.a*sigma.x/sigma.y

  alpha2<-alpha/2
  n2<-n-2
  # critical value of t
  # i.e. upper 100*alpha/2% percentile
  # Pr(t>t.cr)=alpha/2, t~t_{n-2}
  t.cr<-qt(1-alpha2, df=n2)
  part1<-1-pt(t.cr-delta, df=n2)
  part2<-pt(-t.cr-delta, df=n2)

  power<-part1+part2
  res<-list(power=power, delta=delta, s=s, t.cr=t.cr, rho=rho)
  if(verbose)
  { print(res$power) }
  invisible(res)
}

###################################################

tmpPower.SLR<-function(n, power, lambda.a, sigma.x, sigma.y, 
   alpha=0.05, verbose=FALSE)
{
  tmppower<-power.SLR(n=n, lambda.a=lambda.a, sigma.x=sigma.x, 
    sigma.y=sigma.y, alpha=alpha, verbose=verbose)
  res<- tmppower$power-power
  return(res)
}

#####################################################
ss.SLR <- function(power, lambda.a, sigma.x, sigma.y, 
  n.lower = 2.01, n.upper = 1.0e+30, alpha=0.05, verbose=TRUE)
{
  if(n.lower<= 2)
  {
    stop("n.lower must be > 2")
  }
  if(n.lower >= n.upper)
  {
    stop("n.lower must be < n.upper")
  }
  if(power<0 || power > 1)
  {
    stop("power must be in the range [0, 1]")
  }  
  if(lambda.a > sigma.y/sigma.x || lambda.a < -sigma.y/sigma.x)
  {
    stop("lambda.a must be in the range [-sigma.y/sigma.x, sigma.y/sigma.x]")
  }

  res.uniroot<-uniroot(f=tmpPower.SLR, 
      interval=c(n.lower, n.upper),
      power=power, lambda.a=lambda.a, sigma.x=sigma.x, 
      sigma.y=sigma.y, alpha=alpha, verbose=FALSE)

  n.numeric<-res.uniroot$root

  res<-list(n=n.numeric, res.uniroot=res.uniroot)
  if(verbose)
  { print(res$n) }
  invisible(res)
}

############################
tmp.minEffect.SLR<-function(lambda.a, n, power, sigma.x, sigma.y, 
   alpha=0.05, verbose=FALSE)
{

  tmppower<-power.SLR(n=n, lambda.a=lambda.a, sigma.x=sigma.x, 
    sigma.y=sigma.y, alpha=alpha, verbose=verbose)
  res<- tmppower$power-power
  return(res)
}

#####################################################
minEffect.SLR<-function(n, power, sigma.x, sigma.y, 
  alpha=0.05, verbose=TRUE)
{
  if(n <= 2)
  {
    stop("n must be > 2")
  }
  if(power<0 || power > 1)
  {
    stop("power must be in the range [0, 1]")
  }  

  res.uniroot<-uniroot(f=tmp.minEffect.SLR, 
      interval=c(0.0001, sigma.y/sigma.x-0.0001),
      n=n, power=power, sigma.x=sigma.x, 
      sigma.y=sigma.y, alpha=alpha, verbose=FALSE)

  lambda.a<-res.uniroot$root

  res<-list(lambda.a=lambda.a, res.uniroot=res.uniroot)
  if(verbose)
  { print(res$lambda.a) }
  invisible(res)
}

##########################

#########################################################
# for simple linear regression y=theta + lambda * x + epsilon
#########################################################

power.SLR.rho<-function(n, rho2, alpha=0.05, verbose=TRUE)
{
  if(n <= 2)
  {
    stop("n must be > 2")
  }
  if(rho2 > 1 || rho2 <= 0)
  {
    stop("rho2 must be in the range (0, 1]")
  }

  delta=sqrt(n)/sqrt(1/rho2-1)

  alpha2<-alpha/2
  n2<-n-2
  # critical value of t
  # i.e. upper 100*alpha/2% percentile
  # Pr(t>t.cr)=alpha/2, t~t_{n-2}
  t.cr<-qt(1-alpha2, df=n2)
  part1<-1-pt(t.cr-delta, df=n2)
  part2<-pt(-t.cr-delta, df=n2)

  power<-part1+part2
  res<-list(power=power, delta=delta)
  if(verbose)
  { print(res$power) }
  invisible(res)
}

###################################################

tmpPower.SLR.rho<-function(n, power, rho2,
   alpha=0.05, verbose=FALSE)
{
  tmppower<-power.SLR.rho(n=n, rho2=rho2, alpha=alpha, verbose=verbose)
  res<- tmppower$power-power
  return(res)
}

#####################################################
ss.SLR.rho <- function(power, rho2, 
  n.lower = 2.01, n.upper = 1.0e+30, alpha=0.05, verbose=TRUE)
{
  if(n.lower<= 2)
  {
    stop("n.lower must be > 2")
  }
  if(n.lower >= n.upper)
  {
    stop("n.lower must be < n.upper")
  }
  if(power<0 || power > 1)
  {
    stop("power must be in the range [0, 1]")
  }  
  if(rho2 > 1 || rho2 <= 0)
  {
    stop("rho2 must be in the range (0, 1]")
  }

  res.uniroot<-uniroot(f=tmpPower.SLR.rho, 
      interval=c(n.lower, n.upper),
      power=power, rho2=rho2, alpha=alpha, verbose=FALSE)

  n.numeric<-res.uniroot$root

  res<-list(n=n.numeric, res.uniroot=res.uniroot)
  if(verbose)
  { print(res$n) }
  invisible(res)
}


