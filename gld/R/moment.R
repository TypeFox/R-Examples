# Method of Moments estimation for FKML, 
# (c) Sigbert Klinke 2013
# Edits by Robert King (c) 2013,2014,2015
# Licence GPL >= 2
# based on method in 
# Susanna W. M Au-Yeung (2003) Finding Probability Distributions From Moments, Master thesis
# Imperial College of Science, Technology and Medicine (University of London), Department of Computing
# http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.106.6130&rep=rep1&type=pdf
# Additionally, the special cases moments for one of l3,l4 zero by Paul van
# Staden, PhD thesis, University of Pretoria 2013

# To do - integration with fit.fkml
# This function will fit to 4 moments - the fit.fkml version will
# include calculating the moments

# An element of the expression for the moments when l3 and l4 are non-zero
# Thanks to Sigbert Klinke
# k is assumed to be one of 1:4.  Nothing is checked.
vk.nonzero <- function (k, l3, l4) { 
  ret  <- 0
  sign <- 1
  for (j in 0:k) {
    ret  <- ret + sign*choose(k,j)*beta(l3*(k-j)+1, l4*j+1)/(l3^(k-j)*l4^j)
    sign <- -sign
  }
  return (ret)
}

# Calculates the elements of the Moments of the FKML type gld, including
# the limit results for when one of l3 or l4 are zero (or both zero, which
# is the logistic distribution)
# k is assumed to be one of 1:4.  Nothing is checked.
vk <- function(k,l3,l4){
  if (l3==0){
    if(l4==0){ # l3 zero l4 zero - logistic
      return(switch(k,0,(pi^2)/3,0,7*(pi^4)/15))
      } else { # l3 zero l4 non zero
        # psi function is digamma()
        if (k==1){ # See van Staden thesis p 97 
          # k here is 1 there
          # 3 here is k there, therefore -1^k = -1^3 = -1
          # 4 here is j there          
          return((-1*l4)/(l4+1))
        }
        if (k==2){
          element1 = 2*(2*l4^3+l4^2-l4-1)/(l4*(l4+1)*(2*l4+1))
          # Euler's constant, C in van Staden thesis is -digamma(1)
          # psi is digamma().  digamma from base package
          element2 = 2*(-digamma(1)+digamma(l4+2))/(l4*(l4+1))
          return(element1 + element2)
        }
        if (k==3){
         return(m3.1zero.1nonzero(k=3,lj=l4)) # lambda3 is zero, l4 is not
        }
        if (k==4){
          return(m4.1zero.1nonzero(k=3,lj=l4)) # lambda3 is zero, l4 is not
        }
      } 
    } else {
      if(l4==0){# l3 non zero l4 zero
        if (k==1){
          return(l3/(l3+1)) # -1^4
        }
        if (k==2){
          eulC = -digamma(1)
          return(2*(2*l3^3+l3^2-l3-1)/(l3*(l3+1)*(2*l3+1)) +
                   2*(eulC + digamma(l3+2))/(l3*(l3+1)))
        }
        if (k==3){
          return(m3.1zero.1nonzero(k=4,lj=l3)) # l3 is nonzero, l4 is zero
        }
        if (k==4){
          return(m4.1zero.1nonzero(k=4,lj=l3)) # l3 is nonzero, l4 is zero
        }
      } else { # l3 non zero l4 non zero
        return(vk.nonzero(k,l3,l4))
      }
    }
}

m3.1zero.1nonzero <- function(k,lj){
  # This function returns the value of the quantity M3 (from van Staden PhD thesis 2013) in the special cases where one of lambda 3 and lambda 4 is zero and the other is not.  M3 is used in calculating the 3rd moment.  k is the lambda index of the lambda parameter that is zero and lj is the value of the other lambda paramter.
  3*(-1)^k*(12*lj^5+10*lj^4-4*lj^3-lj^2+4*lj+1)/(lj^2*(lj+1)*(2*lj+1)*(3*lj+1)) - 6*
    (-1)^k*(-digamma(1) + digamma(lj+2))/(lj^2*(lj+1)) + 3*
    (-1)^k*(-digamma(1) + digamma(2*lj+2))/(lj^2*(2*lj+1)) + 3*
    (-1)^k*(pi^2/6+(-digamma(1)+digamma(lj+2))^2 - psigamma(lj+2,1))/
    (lj*(lj+1))
}

m4.1zero.1nonzero <- function(k,lj){
  # This function returns the value of the quantity M4 (from van Staden PhD thesis 2013) in the special cases where one of lambda 3 and lambda 4 is zero and the other is not.  M4 is used in calculating the 4th moment.  k is the lambda index of the lambda parameter that is zero and lj is the value of the other lambda paramter.
  zeta3 = 1.20205690315959 # zeta(3) from Rmpfr package, hard-coded in
  # because only zeta(3) is needed 
  eulC = -digamma(1)
  res = 4 * (144*lj^7+156*lj^6-18*lj^5-24*lj^4+7*lj^3-11*lj^2-7*lj-1)/
      (lj^3*(lj+1)*(2*lj+1)*(3*lj+1)*(4*lj+1)) +
    12 * (eulC + digamma(lj+2))/(lj^3*(lj+1)) - 
    12 * (eulC + digamma(2*lj+2))/(lj^3*(2*lj+1)) +
    4  * (eulC + digamma(3*lj+2))/(lj^3*(3*lj+1)) -
    12 * (pi^2/6 + (eulC + digamma(lj+2))^2 - psigamma(lj+2,1))/(lj^2*(lj+1))+
    6  * (pi^2/6 + (eulC + digamma(2*lj+2))^2 - psigamma(2*lj+2,1))/
      (lj^2*(2*lj+1)) + 
    4 * (3 * (pi^2/6 - psigamma(lj+2,1))*(eulC + digamma(lj+2)) +  
      (eulC + digamma(lj+2))^3 + 2*zeta3 + psigamma(lj+2,2) ) / (lj*(lj+1))
  res  
}

fit.fkml.moments.val <- function (moments=c(mean=0, variance=1, 
    skewness=0, kurtosis=3), optim.method="Nelder-Mead",
    optim.control= list(), starting.point = c(0,0)) {
  # optimise.this fitting l3,l4 on the basis of the 3rd & 4th moments 
  # (actually the skewness and kurtosis ratios)
  # par is a vector of length 2, containing lambda3 and lambda4
  # 4 moments exist provided that lambda3 and lambda4 > -0.25
  mean = moments[1]
  variance = moments[2]
  skewness = moments[3]
  kurtosis = moments[4]  
  # Is this the correct thing to do about the the ?skewness missing - does it fix the NaNs?
  optimise.this <- function (par, skewness=moments[3], kurtosis=moments[4]) {
    v1 <- vk(1, par[1], par[2])
    v2 <- vk(2, par[1], par[2])
    v3 <- vk(3, par[1], par[2])
    v4 <- vk(4, par[1], par[2])
    sigmal2sq <- v2-v1^2  # sigmal2sq is sigma^2 hat times lambda2^2, 
    # which is v2 - v1^2.  It is the basis of the denominator of the 
    # skewnesss and kurtosis ratios
    a3 <- (v3-3*v1*v2+2*v1^3)/(sigmal2sq^1.5)
    a4 <- (v4-4*v1*v3+6*v1^2*v2-3*v1^4)/(sigmal2sq^2)
    return((skewness-a3)^2+(kurtosis-a4)^2)
  }
  par       <- starting.point # starting point for l3,l4
  constraint.matrix <- diag(2)
  constraint.vector <- -.25 + .Machine$double.eps
  est       <- constrOptim(theta=par, f=optimise.this, 
      ui=constraint.matrix,ci=constraint.vector,method=optim.method,
      skewness=skewness, kurtosis=kurtosis, control=optim.control)
  lambda    <- c(NA,NA,est$par) # optim has produced estimates for shape parameters
  # Then get lambda2hat, lambda1hat as functions of lambda3hat lambda4hat 
  # and the mean and variance 
  lambda[2] <- sqrt((vk(2, lambda[3], lambda[4])-vk(1, lambda[3], lambda[4])^2)/variance)
  if ((lambda[3]==0)|(lambda[4]==0)|(lambda[3]==lambda[4])){
    # l3 or l4 zero, or l3=l4
    lambda[1] <- mean - (vk(1,lambda[3], lambda[4]))/lambda[2]    
  } else { #l3,l4 nonequal, nonzero case
  lambda[1] <- mean - (vk(1,lambda[3], lambda[4]) - 1/lambda[3] + 1/lambda[4]) /lambda[2]
  }
  names(lambda) <- paste("lambda",1:4,sep="")
  return (lambda)  # estimated lambda parameters 
}

gld.moments <- function(par,type="fkml",ratios=TRUE){
  # par has length 4, the lambda parameters
  # should include an option to only calculate some moments, rather than
  # just insist on 1:4
  if (type != "fkml"){
    if (type != "fmkl") {stop("Only the FKML type is currently implemented")}
  }
  fkml.moments(par=par,ratios=ratios)
}

fkml.moments <- function(par,ratios=TRUE){
  # par is the lambda parameters for the FKML type
  # should include an option to only calculate some moments, rather than
  # just insist on 1:4
  if ((par[4] <= -1)|(par[3] <= -1)) {
    return(c(NA,NA,NA,NA))
  }
  a3 = NA; a4 = NA; m3 = NA; m4 = NA # set to NA in case they aren't calculated
  calc2 = TRUE
  calc3 = TRUE
  calc4 = TRUE
  if ((par[4] <= -1/2)|(par[3] <= -1/2)) {
    calc2 = FALSE
    calc3 = FALSE
    calc4 = FALSE
  }
  if ((par[4] <= -1/3)|(par[3] <= -1/3)) {
    calc3 = FALSE
    calc4 = FALSE
  }
  if ((par[4] <= -1/4)|(par[3] <= -1/4)) {
    calc4 = FALSE
  }
  v1 <- vk(1, par[3], par[4])
  if (calc2){
    v2 <- vk(2, par[3], par[4])
  }
  if (calc3){
    v3 <- vk(3, par[3], par[4])
  }
  if (calc4){
    v4 <- vk(4, par[3], par[4])
  }
  if ((ratios)&(calc3)) { # calc3 implies calc2
    sigmal2sq <- v2-v1^2  # sigmal2sq is sigma^2 hat times lambda2^2, 
    # which is v2 - v1^2.  It is the basis of the denominator of the 
    # skewnesss and kurtosis ratios
    a3 <- (v3-3*v1*v2+2*v1^3)/(sigmal2sq^1.5)
    if (calc4) {a4 <- (v4-4*v1*v3+6*v1^2*v2-3*v1^4)/(sigmal2sq^2)}
  } else {
    if (calc3) {mom3 <- (v3-3*v1*v2+2*v1^3)} # thus ratios
    if (calc4) {mom4 <- (v4-4*v1*v3+6*v1^2*v2-3*v1^4)}
    }
  if (calc2){ sigmasq <- (v2-v1^2)/(par[2]^2) 
    } else {sigmasq <- NA}
  if ((par[3]==0)|(par[4]==0)|(par[3]==par[4])) {
    # l3 or l4 zero, or l3=l4
    mean <- par[1] + v1/par[2]
  } else {
    mean <- par[1] + (v1 -1/par[3]+1/par[4])/par[2]
  }
  if (ratios){
    result <- c(mean,sigmasq,a3,a4)
    names(result) <- c("mean","sigmasq","a3","a4")
  } else {
    result <- c(mean,sigmasq,mom3,mom4)
    names(result) <- c("mean","sigmasq","3rdMom","4thMom")
  }
  result
}
