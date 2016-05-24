predict.tegarch <-
function(object, n.ahead=1,
  initial.values=NULL, n.sim=10000, verbose=FALSE, ...)
{
mFit <- coredata(fitted(object, verbose=TRUE))
n.mFit <- dim(mFit)[1]

if(object$model[1]==1){
  #1-comp model parameters:
  omega <- as.numeric(object$par["omega"])
  phi1 <- as.numeric(object$par["phi1"])
  kappa1 <- as.numeric(object$par["kappa1"])
  if(object$model[2]==0){kappastar<-0}else{kappastar<-as.numeric(object$par["kappastar"])}
  df <- as.numeric(object$par["df"])
  if(object$model[3]==0){skew<-1}else{skew<-as.numeric(object$par["skew"])}

  #initial values:
  if(is.null(initial.values)){
    y.initial <- mFit[n.mFit,"y"]
    lambda.initial <- mFit[n.mFit,"lambda"]
    lambdadagg.initial <- mFit[n.mFit,"lambdadagg"]
  }else{
    y.initial <- initial.values$y
    lambda.initial <- initial.values$lambda
    lambdadagg.initial <- initial.values$lambdadagg
  }
  y2.initial <- y.initial^2
  mueps <- STmean(df, skew=skew)
  dfpluss1 <- df+1
  u.initial <- dfpluss1*(y2.initial+y.initial*mueps*exp(lambda.initial))/(df*exp(2*lambda.initial) * skew^(2*sign( y.initial+mueps*exp(lambda.initial) )) + (y.initial+mueps*exp(lambda.initial))^2) - 1
  gterm.initial <- kappa1*u.initial + kappastar*sign(-y.initial)*(u.initial+1)

  #1-step ahead:
  sigma <- exp(omega + phi1*lambdadagg.initial + gterm.initial)
  stdev <- sigma*sqrt(STvar(df, skew=skew))

  if(n.ahead>1){
    #simulate g-term:
    epsilon <- rST(n.sim, df=df, skew=skew)
    epsilon2 <- epsilon^2
    mueps <- STmean(df, skew=skew)
    eps2muepseps <- epsilon2 - mueps*epsilon
    signeps <- sign(epsilon)
    u <- (df+1)*eps2muepseps/(df*skew^(2*signeps) + epsilon2) - 1
    signnegyupluss1 <- sign(mueps-epsilon)*(u+1)
    gterm.sim <- kappa1*u + kappastar*signnegyupluss1

    for(i in 2:n.ahead){
      Eexpgterm <- rep(NA,i)
      Eexpgterm[1] <- exp(phi1^I(n.ahead-1)*gterm.initial)
      for(j in 2:n.ahead){
        Eexpgterm[j] <- mean(exp(phi1^I(n.ahead-j)*gterm.sim))
      }
      sigma[i] <- exp(omega + phi1^i * lambdadagg.initial) * prod(Eexpgterm)
      stdev[i] <- sigma[i]*sqrt(STvar(df, skew=skew))
    } #end for(i in..)
  } #end if(n.ahead>1)
}else{
  #2-comp model parameters:
  omega <- as.numeric(object$par["omega"])
  phi1 <- as.numeric(object$par["phi1"])
  phi2 <- as.numeric(object$par["phi2"])
  kappa1 <- as.numeric(object$par["kappa1"])
  kappa2 <- as.numeric(object$par["kappa2"])
  kappastar<-as.numeric(object$par["kappastar"])
  df <- as.numeric(object$par["df"])
  if(object$model[3]==0){skew<-1}else{skew<-as.numeric(object$par["skew"])}

  #initial values:
  if(is.null(initial.values)){
    y.initial <- mFit[n.mFit,"y"]
    lambda.initial <- mFit[n.mFit,"lambda"]
    lambda1dagg.initial <- mFit[n.mFit,"lambda1dagg"]
    lambda2dagg.initial <- mFit[n.mFit,"lambda2dagg"]
  }else{
    y.initial <- initial.values$y
    lambda.initial <- initial.values$lambda
    lambda1dagg.initial <- initial.values$lambda1dagg
    lambda2dagg.initial <- initial.values$lambda2dagg
  }
  y2.initial <- y.initial^2
  mueps <- STmean(df, skew=skew)
  dfpluss1 <- df+1
  u.initial <- dfpluss1*(y2.initial+y.initial*mueps*exp(lambda.initial))/(df*exp(2*lambda.initial) * skew^(2*sign( y.initial+mueps*exp(lambda.initial) )) + (y.initial+mueps*exp(lambda.initial))^2) - 1
  gterm1.initial <- kappa1*u.initial
  gterm2.initial <- kappa2*u.initial + kappastar*sign(-y.initial)*(u.initial+1)

  #1-step ahead:
  sigma <- exp(omega+phi1*lambda1dagg.initial+phi2*lambda2dagg.initial+gterm1.initial+gterm2.initial)
  stdev <- sigma*sqrt(STvar(df, skew=skew))

  if(n.ahead>1){
    #simulate g-term:
    epsilon <- rST(n.sim, df=df, skew=skew)
    epsilon2 <- epsilon^2
    mueps <- STmean(df, skew=skew)
    eps2muepseps <- epsilon2 - mueps*epsilon
    signeps <- sign(epsilon)
    u <- (df+1)*eps2muepseps/(df*skew^(2*signeps) + epsilon2) - 1
    signnegyupluss1 <- sign(mueps-epsilon)*(u+1)
    gterm1.sim <- kappa1*u
    gterm2.sim <- kappa2*u + kappastar*signnegyupluss1

    for(i in 2:n.ahead){
      Eexpgterm <- rep(NA,i)
      Eexpgterm[1] <- exp(phi1^I(n.ahead-1)*gterm1.initial + phi2^I(n.ahead-1)*gterm2.initial)
      for(j in 2:n.ahead){
        Eexpgterm[j] <- mean(exp(phi1^I(n.ahead-j)*gterm1.sim + phi2^I(n.ahead-j)*gterm2.sim))
      }
      sigma[i] <- exp(omega + phi1^i * lambda1dagg.initial + phi2^i * lambda2dagg.initial) * prod(Eexpgterm)
      stdev[i] <- sigma[i]*sqrt(STvar(df, skew=skew))
    } #end for(i in..)
  } #end if(n.ahead>1)
} #end if(object$model[1]==...)

#out:
if(verbose){
  out <- as.zoo(cbind(sigma, stdev))
}else{
  out <- as.zoo(stdev)
}
return(out)
}
