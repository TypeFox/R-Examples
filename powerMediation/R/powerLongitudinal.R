########################################
# sample size per group for comparing mean change in longitudinal data with 2 time point
ssLong=function(es, rho=0.5, alpha=0.05, power=0.8)
{
  za=qnorm(1-alpha/2)
  zb=qnorm(power)
  n=4*(1-rho)*(za+zb)^2/es^2

  n.int=ceiling(n)
  return(n.int)
}

# delta = abs(mu1 - mu2)
#  m1 = underlying mean change over time t in group 1
#  m2 = underlying mean change over time t in group 2
# sigma1^2 = variance of baseline values within a treatment group 
# sigma2^2 = variance of follow-up values within a treatment group 
# rho = correlation coefficient between baseline and follow-up values
#   within a treatment group

# full parametrization
ssLongFull=function(delta, sigma1, sigma2, rho=0.5, alpha=0.05, power=0.8)
{
  sigmad2=sigma1^2+sigma2^2-2*rho*sigma1*sigma2
  za=qnorm(1-alpha/2)
  zb=qnorm(power)
  n=2*sigmad2*(za+zb)^2/delta^2

  n.int=ceiling(n)
  return(n.int)
}


# power for comparing mean change in longitudinal data with 2 time point
powerLong=function(es, n, rho=0.5, alpha=0.05)
{
  za=qnorm(1-alpha/2)
  power = pnorm(-za + abs(es)*sqrt(n/(1-rho))/2)

  return(power)
}

# full parametrization
powerLongFull=function(delta, sigma1, sigma2, n, rho=0.5, alpha=0.05)
{
  sigmad2=sigma1^2+sigma2^2-2*rho*sigma1*sigma2
  za=qnorm(1-alpha/2)

  power = pnorm(-za + delta*sqrt(n)/sqrt(2*sigmad2))

  return(power)
}

########################
##########################
##########################
#
# Power calculation for testing if mean changes for 2 groups are the
#     same or not for longitudinal study with more than 2 time points. each subject has n observations
# page 31. Diggle PJ, Liang KY, and Zeger SL (1994). Analysis of Longitundinal Data. Clarendon Press, Oxford


# effect size |beta1A - beta1B|/sigma, where sigma=sd of epsilon_{ij}
  
# m - number of subjects
# n - number of observation per subject
# rho - within subject correlation
  
powerLong.multiTime= function (es, m, nn, sx2, rho = 0.5, alpha = 0.05)
{ 
    sx=sqrt(sx2) 

    #za = qnorm(1 - alpha/2)
    za = qnorm(1 - alpha)
    power = pnorm(-za + abs(es*sx) * sqrt(m*nn/(2*(1 - rho))))
    return(power)
}

ssLong.multiTime= function (es, power, nn, sx2, rho = 0.5, alpha = 0.05)
{ 
    #za = qnorm(1 - alpha/2)
    za = qnorm(1 - alpha)
    zp = qnorm(power)
    m = ceiling(2*(za+zp)^2*(1-rho)/(nn*sx2*es^2))
    return(m)
}

# powerLong.multiTime(es=0.5/sqrt(100), m=196, nn=3, sx2=4.22, rho = 0.5, alpha = 0.05*2)
# ssLong.multiTime(es=0.5/sqrt(100), power=0.8, nn=3, sx2=4.22, rho = 0.5, alpha = 0.05*2)


#######################

# Example 8.33 on page 336 of Rosner (2006)
# n=85
#ssLong(es=5/sqrt(225), rho=0.7, alpha=0.05, power=0.8)
# Example 8.33 on page 336 of Rosner (2006)
# power = 0.8
#powerLong(es=5/sqrt(225), n=85, rho=0.7, alpha=0.05)

# Example 8.34 on page 336 of Rosner (2006)
# power=0.75
#powerLong(es=5/sqrt(225), n=75, rho=0.7, alpha=0.05)


