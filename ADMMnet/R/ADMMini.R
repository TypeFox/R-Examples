

#################################################
#####  Initial values for adaptive methods  #####
#################################################


####################
###  Continuous  ###
####################

IniLm=function(x, y){
  
  N0=nrow(x); p=ncol(x)
  nalambda=10
  
  ### Calculation always based on standardized X and centered y
  tem=scaleC(x)
  xscale=tem$sd; x=tem$x; mx=tem$m
  rm(tem)
  
  my=mean(y); y=y-my
  
  wbeta=rep(1, p)
  ### Lambda path
  lambda_max=maxLambdaLmC(x,y,1.0,wbeta,N0)
  lambda_min=ifelse(N0>=p, lambda_max*0.0001, lambda_max*0.01)
  alambda=lambda_max*(lambda_min/lambda_max)^(c(0:(nalambda-1))/(nalambda-1))
  
  repeat {
    outi=EnetLm(x, y, alpha=0.0, lambda=alambda, keep.beta=T)
    if(!is.null(outi))break
    alambda=alambda*2.0
  }
  
  indexi=ncol(outi$Beta)
  beta0=outi$Beta[, indexi]
  wbeta=1/abs(beta0); sgn=sign(beta0[1:p])
  return(list(wbeta=wbeta, sgn=sgn, lambda=alambda[indexi]))
}


#############
###  Cox  ###
#############

IniCox=function(x, y){
  
  N0=nrow(x);p=ncol(x)
  nalambda=10
  
  prep0=PrepCox(x, y); wbeta=rep(1, p)
  ### Lambda path
  lambda_max=maxLambdaCoxC(prep0$x, prep0$tevent, prep0$N, prep0$nevent, prep0$nevent1, prep0$loc1, prep0$n, 1.0, wbeta, N0)
  lambda_min=ifelse(N0>=p, lambda_max*0.0001, lambda_max*0.01)
  alambda=lambda_max*(lambda_min/lambda_max)^(c(0:(nalambda-1))/(nalambda-1))
  
  repeat {
    outi=EnetCox(x, y, alpha=0.0, lambda=alambda, keep.beta=T)
    if(!is.null(outi))break
    alambda=alambda*2.0
  }
  
  indexi=ncol(outi$Beta)
  beta0=outi$Beta[, indexi]
  wbeta=1/abs(beta0); sgn=sign(beta0[1:p])
  return(list(wbeta=wbeta, sgn=sgn, lambda=alambda[indexi]))
}




