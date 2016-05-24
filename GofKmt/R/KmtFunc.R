#' Performs goodness-of-fit test through Khmaladze matringale transformation
#'@param X  - Random sample of n observations
#'@param F0 - Name of null distribution: Normal, Cauchy, or Logistic 
#'@return TestStat - Test statistic obtained through Khmaladze martingale transformation 
#'@return CritValue - Vector of critical values for the level of 0.01, 0.025, 0.05, and 0.10 
#'@return Mu - Maximum likelihood estimator of location parameter mu 
#'@return Sigma - Maximum likelihood estimator of scale parameter sigma
#'@examples
#'####################
#'n = 10
#'Sample = rnorm(n, 1,3)    # Generate a random sample of n observations from N(1,3)
#'KMT_Result = KhmaladzeTrans(Sample, "Normal")
#'KMT_TestStat = KMT_Result$TestStat
#'KMT_CriticalValue = KMT_Result$CritValue
#'KMT_Muhat = KMT_Result$Mu
#'KMT_Sigmahat = KMT_Result$Sigma
#'#####################
#'
#'#####################
#'n = 10
#'Sample = rcauchy(n, 0,2)  # Generate a random sample of n observations from Cauchy distribution
#'KMT_Result = KhmaladzeTrans(Sample, "Cauchy")
#'KMT_TestStat = KMT_Result$TestStat
#'KMT_CriticalValue = KMT_Result$CritValue
#'KMT_Muhat = KMT_Result$Mu
#'KMT_Sigmahat = KMT_Result$Sigma
#'#####################



#'@references
#'[1] E.V. Khmaladze, H.L. Koul (2004). Martingale Transforms Goodness-of-fit tests in regression models. Ann. Statist., 32. 995-1034 
#'@references
#'[2] E.V. Khmaladze, H.L. Koul (2009). Goodness-of-fit problem for errors in nonparametric regression: Distribution free approach. Ann. Statist., 37(6A) 3165-3185. 
#'@importFrom "stats" "optim" "dnorm" "integrate" "median" "optimize" "pnorm" "qcauchy" "qlogis" "qnorm" "quantile" "sd" 
#'@importFrom "Rsolnp" "solnp"
#'@export

KhmaladzeTrans = function(X, F0){
  
  SampleSize = length(X)
  
  if(F0 == "Normal"){
    muhat = mean(X)
    sigmahat = sd(X)
    
  }else if(F0 == "Cauchy"){
    
    x01 = median(X)
    
    quan3 = quantile(X, 0.75)
    quan1 = quantile(X, 0.25)
    
    x02 = (quan3[[1]]-quan1[[1]])/2

    ####################
    x0 = c(x01, x02)
    LBV = c(-Inf, 0)
    UBV = c(Inf, Inf)

    MRe = solnp( pars=x0, fun=MLECauchyObjFunction(X), eqfun=NULL, eqB=NULL, ineqfun=NULL, ineqLB = NULL, 
                 ineqUB = NULL, LB = LBV, UB = UBV)
    
    CauchyMle = MRe$pars
    
    muhat = CauchyMle[1]
    sigmahat = CauchyMle[2]
    
  }else if(F0 == "Logistic"){
    
    x01 = median(X)
    x02 = sqrt(3)/pi*sd(X)
    
    startVal = c(x01, x02)
    mle_method = optim(startVal, mleLogis(X))
    mle_par = mle_method$par
    
    muhat = mle_par[1]
    sigmahat = mle_par[2]
    
  }else{
    message("Name of null distribution is not valid.")
    stop()
  }
  
  StandX = rep(0, times=SampleSize) 
  
  for(i in 1:SampleSize){
    StandX[i] = (X[i]-muhat)/sigmahat
    
  }
  
  d_alpha= c(2.80705, 2.49771, 2.24241, 1.959, 1.64)
  
  if(F0 == "Normal"){
    KhmResult = SupWNormal(StandX)
  }else if(F0 == "Cauchy"){
    KhmResult = SupWCauchy(StandX)
  }else if(F0 == "Logistic"){
    KhmResult = SupWLogistic(StandX)
  }  
  
  lst = list(TestStat=KhmResult$objective, CritValue = c("1%: 2.80705", "2.5%: 2.49771", "5%: 2.24241", "10%: 1.959"), 
            Mu = muhat, Sigma = sigmahat )
  return(lst)
}



################### Normal Part

fn=function(x){
  dnorm(x,0,1)
}

Fn = function(x){
  pnorm(x,0,1)
}

Fn1 = function(x){
  1-pnorm(x,0,1)
}


cx = function(x){
  -(x^2+1)*(fn(x))^2 + (x^3+3*x)*fn(x)*Fn1(x) + 2*(Fn1(x))^2
}

dx = function(x){
  2*(Fn1(x))^3 + (x^3+3*x)*fn(x)*( Fn1(x))^2 - (2*x^2+3)*Fn1(x)*( fn(x))^2 + x*( fn(x))^3
}

c22=function(x){
  1/cx(x)
}

d22=function(x){
  1/dx(x)
}


D1=function(x){
  2*fn(x)* ( (Fn1(x))^2 + x*fn(x)*Fn1(x) - (fn(x))^2 )/dx(x)
}

D2 = function(x){
  ans = 4*x* ( Fn1(x) )^4 +(2*x^4+8*x^2 -2)*fn(x)*( Fn1(x) )^3 + (x^5-7*x)*( fn(x) )^2*( Fn1(x) )^2 -
    (2*x^4 + 3*x^2 - 1)*( fn(x) )^3*Fn1(x) + (x^3+x)*( fn(x) )^4
  return(ans*fn(x)*c22(x)*d22(x))   
}

D3 = function(x){
  
  ans = (2*x^2-2)* ( Fn1(x) )^4 +( x^5 + 2*x^3 -9*x)*fn(x)*( Fn1(x) )^3 - (4*x^4 + 9*x^2 -5)*( fn(x) )^2*( Fn1(x) )^2 +
    (5*x^3 + 9*x )*( fn(x) )^3*Fn1(x) - (2*x^2 + 2)*( fn(x) )^4
  
  return(ans*fn(x)*c22(x)*d22(x))   
}


ItgofIntNormal = function(ri){
  Dual = function(y){
    ans = D1(y) + ri*D2(y)+(ri^2-1)*D3(y)
    return(ans)
  }
  return(Dual)
}


IntNormal = function(ri,t){
  upValVec = c(ri, qnorm(t, 0,1))
  upVal = min(upValVec)
  ans = integrate( ItgofIntNormal(ri), lower=-Inf, upper = upVal)$value
  return(ans)
  
}


WNormal = function(r){
  
  InsideW1n = function(t){
    x = qnorm(t, 0,1)
    n=length(r)
    tempsum = 0
    for (i in 1:n){
      
      tempsum = tempsum + ( (r[i]<=x) - IntNormal(r[i],t)  )
      
    }
    tempsum = tempsum/sqrt(n)
    return(abs(tempsum))	
  }
  return(InsideW1n)
  
}

SupWNormal = function(r){
  
  tmin = optimize(WNormal(r),  lower=0, upper =1, maximum=TRUE )
  
}


##############################################



#################################### Cauchy

##################### Cauchy MLE

MLECauchyObjFunction = function(XVec){
  nLength = length(XVec)
  
  Dual = function(x){
    tempsum = -nLength*log(pi*x[2])
    for(i in 1:nLength){
      tempsum = tempsum - log(1+( (XVec[i]-x[1] ) / x[2]  )^2)
    }
    return(-tempsum)
  }
  return(Dual)
}


###################  
dc = function(x){
  (pi-2*atan(x))
}

bc = function(x){
  (1+x^2)*dc(x)-2*x
  
}

ac = function(x){
  (1+x^2)*dc(x)^2 + 2*x*dc(x)-8
}



KhmIntValCauchy = function(ri,t){
  
  upValVec = c(ri, qcauchy(t, 0,1))
  upVal = min(upValVec)
  
  funcVal = integrate( ItgofIntCauchy(ri), lower=-Inf, upper = upVal)$value
  return(funcVal)
}



ItgofIntCauchy = function(ri){
  Dual = function(x){
    ##################
    
    num = bc(x)*( 2*bc(x) -(16*ri+8*(ri^2-1)*x )/(1+ri^2)  ) +
      (16*ri*x)/(1+ri^2)*( (1+x^2)*dc(x)^2-4 ) +
      4*(ri^2-1)/(1+ri^2)*( (x^4-1)*dc(x)^2 - 2*x*(1+x^2)*dc(x)+8    )
    
    
    den = (1+x^2)*ac(x)*bc(x)
    
    
    dualVal = num/den
    
    ########    
    
    #dualVal = I1(x)-I2(x)-I3(ri,x)+I4(ri,x)+I5(ri,x)
    
    return(dualVal)
  }
  return(Dual)
}
#####################################



WCauchy = function(r){
  
  InsideW1n = function(t){
    x = qcauchy(t, 0,1)
    n=length(r)
    tempsum = 0
    for (i in 1:n){
      
      tempsum = tempsum + ( (r[i]<=x) - KhmIntValCauchy(r[i],t)  )
      
    }
    tempsum = tempsum/sqrt(n)
    return(abs(tempsum))	
  }
  return(InsideW1n)
  
}

SupWCauchy = function(r){
  
  tmin = optimize(WCauchy(r),  lower=0.01, upper =1, maximum=TRUE )
  
}



########################### Logistic

#######################

mleLogis = function(X){
  
  nLeng = length(X)
  Dual = function(param){
    tempsum=0
    for (i in 1:nLeng){
      normx = (X[i]-param[1])/param[2]
      tempsum = tempsum + normx - log(param[2]) -2*log(1+exp(normx))
    }
    return(-tempsum)
  }
  return(Dual)
  
}
###########


######################

fl = function(x){
  exp(x)/( 1+exp(x) )^2
}


Fl = function(x){
  exp(x)/(1+exp(x))
}


g1 = function(x){
  log(1+exp(x))/3 - fl(x)* (x*(3+exp(2*x))+(1+exp(x)) )/(3*(1+exp(x)))
}


Weird = function(x){
  x^2*((1-exp(x))/(1+exp(x)))^2*fl(x)
  
}

Gam333 = function(x){
  
  ans = integrate( Weird, lower=x, upper = 600)$value 
  return(ans) 
}

g2 = function(x){
  Fl(x)-2*x*fl(x)+Gam333(x)-1  
}




k1 = function(x){
  3*exp(x)*(1+exp(x))^3*g2(x) - 3*x*exp(x)*(1+exp(x))^3*g1(x)
  
}

k2 = function(x){
  -3*exp(x)*(1+exp(x))^3*g1(x) + x*exp(x)*(1+3*exp(2 *x))
  
}


cofx = function(x){
  g1 = g1(x)
  g2 = g2(x)
  funcVal = 1/( 3*(1+3*exp(2*x))*(1+exp(x))^3*g2 - 9*(1+exp(x))^6*g1^2 )
  return(funcVal)
  
}

A12A22A21=function(x){
  return(9*cofx(x)*exp(2*x)*(1+exp(x))^2* 
           ( 3*(1+exp(x))^3*g2(x)- 6*x*(1+exp(x))^3*g1(x) +x^2*(1+3*exp(2*x))) )
}

B11 = function(x){
  funcVal = 1/(  3*(1+exp(x))^2  - A12A22A21(x))
  return(funcVal)
}


KhmIntValLogistic = function(ri,t){
  
  upValVec = c(ri, qlogis(t, 0,1))
  upVal = min(upValVec)
  
  funcVal = integrate( ItgofIntLogistic(ri), lower=-Inf, upper = upVal)$value
  
  return(funcVal)
}


ItgofIntLogistic = function(ri){
  Dual = function(x){
    
    dualVal = 3*exp(x)*(1+exp(x))*B11(x)*(1+3*cofx(x)*( (1-exp(x))*k1(x) + (1+exp(x)+x-x*exp(x) )*k2(x) )  )+ 
      ( -(1-exp(ri))/(1+exp(ri))  ) * ( (  -9*exp(x)*(1+exp(x))^2* B11(x)*cofx(x) * k1(x) ) -  
                                          3*exp(x)*cofx(x)*(3*(1+exp(x))^3*(1-exp(x) )*g2(x)- 3*(1+exp(x))^3*g1(x)* (1+exp(x)+x-x*exp(x)) ) -
                                          B11(x)*27*exp(x)* (1+exp(x))^2* (cofx(x))^2 * ( (1-exp(x))*(k1(x))^2 +(1+exp(x)+x-x*exp(x))*k1(x)*k2(x)   )  )+  
      ( -1 - ri*(1-exp(ri))/(1+exp(ri))  ) * (  (  -9*exp(x)*(1+exp(x))^2* B11(x)*cofx(x) * k2(x) ) -  
                                                  3*exp(x)*cofx(x)*( -3*(1+exp(x))^3*(1-exp(x) )*g1(x) + (1+3*exp(2*x))* (1+exp(x)+x-x*exp(x)) )- 
                                                  B11(x)*27*exp(x)* (1+exp(x))^2* (cofx(x))^2 * ( (1-exp(x))*k1(x)*k2(x) +(1+exp(x)+x-x*exp(x))* (k2(x))^2   ) )    
    
    return(dualVal)
  }
  return(Dual)
}
#####################################



WLogistic = function(r){
  
  InsideW1n = function(t){
    x = qlogis(t, 0,1)
    n=length(r)
    tempsum = 0
    for (i in 1:n){
      tempsum = tempsum + ( (r[i]<=x) - KhmIntValLogistic(r[i],t)  )
    }
    tempsum = tempsum/sqrt(n)
    return(abs(tempsum))	
  }
  return(InsideW1n)
  
}

SupWLogistic = function(r){
  
  tmin = optimize(WLogistic(r),  lower=0, upper=1, maximum=TRUE )
  
}










