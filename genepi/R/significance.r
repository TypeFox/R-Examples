# R code to estimate two-stage design parameters
#
# revised: march 5, 2003
#
# jaya m. satagopan (satagopj@mskcc.org)
#
# reference: Satagopan JM, Elston RC (2003): Optimal two-stage genotyping in population-based
#                 association studies. Genetic Epidemiology 25: 149-157.
#
#
# how to use these functions:
#
# First read this file in R using the "source" command.
# For example, if you name this file as "twostage.s", then
# in R, use the command - source("twostage.s"). Remember
# to give the appropriate path. You may want to say -
# source("c:/twostagecode/twostage.s") - if you have saved
# this file under the directory c:\twostagecode. 
#
# The main function is "find.maximum.power". Run this 
# function to get the required results. 
#
# Example to run this function:
# Suppose we have m=25, D=1, alpha (overall significance level)=0.05, power (of onestage method)=0.80, 
# cost fraction=0.75, and alpha1 (significance level at Stage 1) = 0.20. Then, use the following command:
# 
#       find.maximum.power(m=25, Dvalue=1, alpha=0.05, power=0.80, cost.frac=0.75, alpha1=0.20)
#
# Note: The code given here does not perform the grid search to find the optimum power and optimum
# cost. This must be coded by the user. It can be done easily by using the function find.maximum.power
# iteratively. In fact, if you give a vector of values for alpha1, the function find.maximum.power will output
# the values of alpha1, beta1, and alpha2 at which the maximum power occurs (and the value of the
# maximum power), for a given cost. So, the iteration is only required to iterate over various values of
# cost.frac. 
#
# Should you have any questions regarding these functions, please contact 
# Jaya M. Satagopan at 646-735-8122 or satagopj@mskcc.org
#

###############################
# main function to generate the required results
###############################
find.maximum.power <- function(m=100, Dvalue=1, alpha=0.05, power=0.80, cost.frac=0.25, alpha1 =0.25){
  vec.len <- length(alpha1)
  beta1 <- rep(0,vec.len)
  alpha2 <- rep(0,vec.len)
  power.twostage <- rep(0, vec.len)
  power.twostage.exact <- rep(0, vec.len)
  for(i in 1:vec.len){
    temp <- getestimates(m=m, Dvalue=Dvalue, alpha=alpha, power=power, cost.frac=cost.frac, alpha1.start=alpha1[i])
    beta1[i] <- temp$beta1
    alpha2[i] <- temp$alpha2
    power.twostage[i] <- temp$power.twostage
  }
  max.loc <- (1:vec.len)[power.twostage >= max(power.twostage)]
    return(list(m=m, D=Dvalue, alpha=alpha, power=power, cost.frac=cost.frac,
                 opt.alpha1=alpha1[max.loc], opt.pow1=1-beta1[max.loc], opt.alpha2=alpha2[max.loc], opt.power=power.twostage[max.loc]))
}


####################################
# function to calculate the estimates
####################################
getestimates <- function(m=100, Dvalue=1, alpha=0.05, power=0.80, cost.frac=0.55, 
                                        alpha1.start=0.10, beta1.start=0.000001, alpha2.start=0.00000001){

    alpha1 <- alpha1.start

    beta1 <- get.beta1(m=m, Dvalue=Dvalue, alpha=alpha, power=power, alpha1=alpha1.start, 
                                   cost.frac=cost.frac,beta1.start=beta1.start)
    alpha2 <- get.alpha2(m=m, alpha=alpha, power=power, 
                                      alpha1=alpha1, beta1=beta1, alpha2.start=alpha2.start)

    twostage.power.exact <- twostage.power.exact(m=m, alpha=alpha, power=power, alpha1=alpha1, beta1=beta1, alpha2=alpha2)

    return( list(m=m, alpha=alpha, power=power,
                 cost.frac=cost.frac, Dvalue=Dvalue, 
                 alpha1=alpha1, beta1=beta1,
                 power1=1-beta1, alpha2=alpha2,
                 power.twostage=twostage.power.exact))

}

############################################
# estimate beta1 using newton raphson method
############################################
get.beta1 <- function(m=100, Dvalue=1, alpha=0.05, power=0.80, alpha1=0.001, 
                                   cost.frac=0.8, beta1.start=0.000001, threshold=0.000001){
  epsilon <- 1
  iteration <- 0
  beta1.start <- beta1.start
  while(epsilon > threshold){
      f.beta1 <- ( m*(qnorm(1-alpha1)+qnorm(1-beta1.start))^2 + 
                     ( (m-Dvalue)*alpha1 + Dvalue*(1-beta1.start) ) * 
                       ( (qnorm(1-alpha/m)+qnorm(power))^2 - (qnorm(1-alpha1)+qnorm(1-beta1.start))^2 ) ) /
                     ( m * (qnorm(1-alpha/m)+qnorm(power))^2 ) - cost.frac
      f.prime.beta1 <- 2*(qnorm(1-alpha1)+qnorm(1-beta1.start))/dnorm(qnorm(1-beta1.start)) * 
                              ( (m-Dvalue)*alpha1 + Dvalue*(1-beta1.start) - m) -
                              Dvalue*( (qnorm(1-alpha/m)+qnorm(power))^2 - (qnorm(1-alpha1)+qnorm(1-beta1.start))^2 )
      f.prime.beta1 <- f.prime.beta1 / ( m * (qnorm(1-alpha/m)+qnorm(power))^2 ) 
      beta1.new <- beta1.start - f.beta1/f.prime.beta1
      if(beta1.new < 0){
        print("negative beta1 estimate. convergence problem!! change the starting value of beta1")
        break
      }
      iteration <- iteration + 1
      if(iteration > 20){
        print("convergence problem! change the starting value of beta1")
        break
      }
      epsilon <- abs(beta1.new - beta1.start)
      beta1.start <- beta1.new
  }
  return(beta1.new)
}

#############################################
# estimate alpha2 using newton raphson method
#############################################
get.alpha2 <- function(m=100, alpha=0.05, power=0.80, 
                                    alpha1=0.001, beta1=0.03, alpha2.start=0.000000001, threshold=0.000001){
  epsilon <- 1
  iteration <- 1
  alpha2.start <- alpha2.start

  while(epsilon > threshold){
    f.alpha2 <- integrate(integrand, lower=qnorm(1-alpha1), upper=10, m=m, alpha=alpha, power=power, alpha1=alpha1, beta1=beta1, alpha2=alpha2.start)$value - alpha/m
    f.prime.alpha2 <- integrate(integrand.prime, lower=qnorm(1-alpha1), upper=10, m=m, alpha=alpha, power=power, alpha1=alpha1, beta1=beta1, 
                                           alpha2=alpha2.start)$value
    alpha2.new <- alpha2.start - f.alpha2/f.prime.alpha2
    epsilon <- abs(alpha2.new - alpha2.start)
    alpha2.start <- alpha2.new
  }
  return(alpha2.new)
}

#######################################
# integrand for calculating alpha2
#######################################
integrand <- function(z.value, m, alpha, power, alpha1, beta1, alpha2){
  temp1 <- qnorm(1-alpha/m) + qnorm(power)
  temp2 <- qnorm(1-alpha1) + qnorm(1-beta1)
  numerator <- qnorm(1-alpha2) * temp1 -z.value*temp2
  denominator <- sqrt( temp1^2 - temp2^2 )
  (1 - pnorm(numerator/denominator) ) *dnorm(z.value)
}

#########################################################
# integrand for the derivative of alpha/m with respect to alpha2
#########################################################
integrand.prime <- function(z.value, m, alpha, power, alpha1, beta1, alpha2){
  temp1 <- qnorm(1-alpha/m) + qnorm(power)
  temp2 <- qnorm(1-alpha1) + qnorm(1-beta1)
  numerator <- qnorm(1-alpha2) * temp1 - z.value * temp2
  denominator <- sqrt(temp1^2-temp2^2)
  1/dnorm(qnorm(1-alpha2)) * temp1/sqrt(temp1^2-temp2^2) * dnorm(numerator/denominator) * dnorm(z.value)
}

##################################
# power P of the two-stage design
##################################
twostage.power.exact <- function(m=100, alpha=0.05, power=0.80, alpha1=0.10, beta1=0.05, alpha2=0.0005){
  power.value.exact <- integrate(integrand.power.exact, lower=qnorm(1-alpha1), upper=10, m=m, alpha=alpha, power=power, alpha1=alpha1, beta1=beta1, alpha2=alpha2)$value
  return(power.value.exact)
}

##################################
# integrand of the power P (see above function twostage.power.exact)
##################################
integrand.power.exact <- function(z.value, m, alpha, power, alpha1, beta1, alpha2){
  temp1 <- qnorm(1-alpha/m) + qnorm(power)
  temp2 <- qnorm(1-alpha1) + qnorm(1-beta1)
  temp3 <- temp1^2 - temp2^2
  numerator <- qnorm(1-alpha2) * temp1 - z.value*temp2 - temp3
  denominator <- sqrt(temp3)
  integrand.value <- ( 1 - pnorm(numerator/denominator) ) * dnorm(z.value - temp2)
  return(integrand.value)
}


