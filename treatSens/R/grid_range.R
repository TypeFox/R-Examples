#############
#Find boundaries for grid search
#############

###
#Divide and conquer algorithm
#Finds where tau=0 is crossed on diagonal
#arguments:
#x1, x2: Z endpoints of interval
#y1, y2: Y endpoints of interval
#tau1,tau2: treatment effect estimates at endpoints
#fn.call: sensitivity analysis function call (will call one of fit_GLM_sens or fit_BART_sens)
###

DandCsearch <- function(x1, x2, y1, y2, tau1, tau2, fn.call) {
  #change coefficients in function call
  fn.call[8] = (y1+y2)/2
  fn.call[9] = (x1+x2)/2
  
  aaa = eval(fn.call)
  tauM = aaa$sens.coef
  
  if(abs(tauM) < aaa$sens.se/3){
    return(list(rZ = (x1+x2)/2, rY = (y1+y2)/2, tau = tauM))
  }else{
    if(sign(tau1)==sign(tauM)) 
      DandCsearch((x1+x2)/2, x2, (y1+y2)/2, y2, tauM, tau2, fn.call)
    else
      DandCsearch((x1+x2)/2, x1, (y1+y2)/2, y1, tauM, tau1, fn.call)
  }
}

###
#Find range for grid search
#arguments:
#extreme.cors: matrix of extreme possible correlations
#Xcoef.plot: plotting coordinates for X values
#Y: response
#Z: treatment
#X: matrix of covariates
#Y.res: residuals of Y from null fit (Y~Z+X)
#Z.res: residuals of Z from null fit (Z~X)
#resp.family, trt.family, U.model: see GLM.sens or BART.sens
#sgnTau0: the sign of the estimated treatment effect in the null model
###

grid.search <- function(extreme.cors, zero.loc, Xcoef.plot, Y, Z, X, 
                        Y.res, Z.res, v_Y, v_Z, theta, sgnTau0, control.fit, sensParam) {
  if(!is.null(control.fit$resp.family)){
    fname <- "fit.treatSens"
  }else{
    fname <- ifelse(is.null(control.fit$g), "fit.treatSens.BART", "fit.LMER.sens")
  }
  fn.call <- call(fname, sensParam = sensParam, Y=Y, Z=Z, Y.res=Y.res, Z.res=Z.res, X=X, zetaY=NA, zetaZ=NA, v_Y = v_Y, v_Z = v_Z,
                  control.fit = control.fit, theta = theta)
  
  tau = rep(NA, 3)
  #  calculate new tau for (extreme Y)*(0, midpoint, extreme Z) to see if 0 is crossed (and in which half)
  #		note that extreme depends on sign of tau: ++ for positive tau, +- for negative
  if(sgnTau0 == -1) {
    rY = c(0,
           extreme.cors[1,2]/2, #midpoint
           extreme.cors[1,2]*.95) #extreme Y
    rZ = c(0,
           extreme.cors[2,1]/2, #midpoint
           extreme.cors[2,1]*.95) #extreme Z
    
  }else{
    rY = c(0,
           extreme.cors[1,2]/2, #midpoint
           extreme.cors[1,2]*.95) #extreme Y
    
    rZ = c(0,
           extreme.cors[2,2]/2, #midpoint
           extreme.cors[2,2]*.95) #extreme Z
  }
  for(i in 1:3) {
    fn.call[8] = rY[i]
    fn.call[9] = rZ[i]
    aaa = eval(fn.call)
    tau[i] =aaa$sens.coef	
  }
  
  if(sign(tau[1]) == sign(tau[3])) {
    Ymax = rY[3]
    Zmax = rZ[3]
  }else {
    loc0 <- DandCsearch(rZ[2], rZ[sign(tau)!=sign(tau[2])], rY[2], rY[sign(tau)!=sign(tau[2])], tau[2], tau[sign(tau)!=sign(tau[2])], fn.call) 
    Zmax = sign(rZ[3])*min(abs(rZ[3]),abs(1/zero.loc*loc0$rZ))
    
    #update Zmax so that it include all plots.
    if(!is.null(X)){
      if ((sgnTau0 == 1) && (all(Xcoef.plot[,1]<extreme.cors[2,2]*.95))) {
        Zmax = ifelse(Zmax>max(Xcoef.plot[,1]),Zmax,max(Xcoef.plot[,1]))}
      if ((sgnTau0 == -1) && (all(Xcoef.plot[,1]>extreme.cors[2,1]*.95))) {
        Zmax = ifelse(Zmax<max(Xcoef.plot[,1]),Zmax,max(Xcoef.plot[,1]))}
    }
    Ymax = min(rY[3],sqrt(v_Y)/(1-Zmax^2/v_Z))
    if(control.fit[3] == "binomial") Ymax = rY[3]
  }
  
  ######Do it again for opposite side
  if(sgnTau0 == 1) {
    rY = c(0,
           extreme.cors[1,2]/2, #midpoint
           extreme.cors[1,2]*.95) #extreme Y
    rZ = c(0,
           extreme.cors[2,1]/2, #midpoint
           extreme.cors[2,1]*.95) #extreme Z
    
  }else{
    rY = c(0,
           extreme.cors[1,2]/2, #midpoint
           extreme.cors[1,2]*.95) #extreme Y
    
    rZ = c(0,
           extreme.cors[2,2]/2, #midpoint
           extreme.cors[2,2]*.95) #extreme Z
  }
  for(i in 1:3) {
    fn.call[8] = rY[i]
    fn.call[9] = rZ[i]
    aaa = eval(fn.call)
    tau[i] =aaa$sens.coef
    
  }
  
  if(sign(tau[1]) == sign(tau[3])) {
    Ymin = rY[3]
    Zmin = rZ[3]
  }else {
    loc0 <- DandCsearch(rZ[2], rZ[sign(tau)!=sign(tau[2])], rY[2], rY[sign(tau)!=sign(tau[2])], tau[2], tau[sign(tau)!=sign(tau[2])], fn.call) 
    Zmin = sign(rZ[3])*min(abs(rZ[3]),abs(1/zero.loc*loc0$rZ))
    
    #update Zmin so that it include all plots.
    if(!is.null(X)){
      if ((sgnTau0 == 1) && all(Xcoef.plot[,1]>extreme.cors[2,1]*.95)) {
        Zmin = ifelse(Zmin<min(Xcoef.plot[,1]),Zmin,min(Xcoef.plot[,1]))}
      if ((sgnTau0 == -1) && all(Xcoef.plot[,1]<extreme.cors[2,2]*.95)) {
        Zmin = ifelse(Zmin>min(Xcoef.plot[,1]),Zmin,min(Xcoef.plot[,1]))}
    }
    Ymin = max(rY[3],sqrt(v_Y)/(1-Zmin^2/v_Z))
    if(control.fit[3] == "binomial") Ymax = rY[3]
  }
  
  Z.range = c(min(Zmin, Zmax), max(Zmin, Zmax))
  Y.range = c(0, min(Ymin, Ymax))
  
  return(ranges = rbind(Y.range, Z.range))
}
