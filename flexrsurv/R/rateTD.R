rateTD_beta0alphabeta<- function(T, iT, gamma0, Zbeta0, Zalphabeta, 
                      Spline_t0=SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t0=TRUE,
                      Spline_t =SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t=TRUE, ...){
  # compute the contribution of the time dependent variables to the rate 
  # of relative survival model for patient iT with Zalphabeta[iT, ]
  # at a vector of T (useful to compute numerical integration 
  # spline bases for baseline hazard
  # all T basis for the NPH/td effects are the same (Spline_t)
  # Zbeta0     = X %*% beta0 
  # Zalphabeta = f(Z,alpha) %*% beta 

  #  YT0 <- bs(T, knots=Knots_t0, intercept=Intercept_t0, degree=degree_t0, Boundary.knots =  Boundary.knots_t0)
    YT0 <- fevaluate(Spline_t0, T, intercept=Intercept_t0, outer.ok=TRUE)
  # spline bases for each TD effect
#  YT  <- bs(T, knots=Knots_t, intercept=Intercept_t, degree=degree_t, Boundary.knots =  Boundary.knots_t)
    if(!is.null(Zbeta0)) {
      if(!is.null(Zalphabeta)) {
        exp(YT0 %*% gamma0 +
            fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE) %*% Zbeta0[iT,] +
            fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE) %*% Zalphabeta[iT,])
      }
      else {
        exp(YT0 %*% gamma0 +
            fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE) %*% Zbeta0[iT,] )
      }
    }
    else if(!is.null(Zalphabeta)) {
        exp(YT0 %*% gamma0 +
            fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE) %*% Zalphabeta[iT,])
      }
}

rateTD_gamma0alphabeta<- function(T, iT, gamma0, Zalphabeta, 
                      Spline_t0=SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t0=TRUE,
                      Spline_t =SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t=TRUE, ...){
  # compute the contribution of the time dependent variables to the rate 
  # of relative survival model for patient iT with Zalphabeta[iT, ]
  # at a vector of T (useful to compute numerical integration 
   # all T basis for the NPH/td effects are the same (Spline_t)
  # Zalphabeta = X %*% beta0 + f(Z,alpha) %*% beta 

  # spline bases for baseline hazard
    YT0 <- fevaluate(Spline_t0, T, intercept=Intercept_t0, outer.ok=TRUE)
  # spline bases for each TD effect
    YT <- fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE)

  # returned value 
  exp(YT0 %*% gamma0 + YT %*% Zalphabeta[iT,])
  
}

rateTD_alphabeta<- function(T, iT, Zalphabeta, 
                      Spline_t =SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t=TRUE, ...){
  # compute the contribution of the time dependent variables to the rate 
  # of relative survival model for patient iT with Zalphabeta[iT, ]
  # at a vector of T (useful to compute numerical integration 
  # all T basis for the NPH/td effects are the same (Spline_t)
  # Zalphabeta = X %*% beta0 + f(Z,alpha) %*% beta 
 
  # spline bases for each TD effect
    YT <- fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE)

  # returned value 
  exp(YT %*% Zalphabeta[iT,])
  
}

rateTD_gamma0<- function(T, iT, gamma0, 
                      Spline_t0=SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t0=TRUE, ...){
  # compute the contribution of the baseline hazard rate to the rate 
  # of relative survival model for patient iT with Zalphabeta[iT, ]
  # at a vector of T (useful to compute numerical integration 
 
  # spline bases for baseline hazard
    YT0 <- fevaluate(Spline_t0, T, intercept=Intercept_t0, outer.ok=TRUE)
  # returned value

  exp(YT0 %*% gamma0)
  
}


rateTD<- function(T, iT, ... ){
  # compute the contribution of the baseline hazard rate when there are no time dependent effect
  # and when the baseline hazard is null
 
  # returned value

  rep(1, length(T))
  
}



# computes the contribution of time dependent termes in the rate (baseline, NPH, NPHNLL and WCE effects 
rateTD_gamma0alphabetaeta0<- function(T, iT,
                                      fromT, toT, fail, FirstId,
                                      gamma0, Zalphabeta,
                                      nW, W, eta0, iWbeg, iWend,
                      Spline_t0=SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t0=TRUE,
                      Spline_t =SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t=TRUE,
                                      ISpline_W, Intercept_W=rep(TRUE,nW), ...){
  # compute the contribution of the time dependent variables to the rate 
  # of relative survival model for line iT with Zalphabeta[iT, ]
  # at a vector of T (useful to compute numerical integration
  # fromT : begining of the time intervals
  # toT   : end of the time intervals
  # fail  : fail = 1 if event, 0 if censored
  # all T basis for the NPH/td effects are the same (Spline_t)
  # Zalphabeta = X %*% beta0 + f(Z,alpha) %*% beta 
  # FirstId : all lines in FirstId[iT]:iT of fromT, toT, fail and Zalphabeta comes from the same individual 
  # nW number of cols in W (number of WCE effects
  # W matrix of exposure INCREMENT variables W[FirstId[iT]:iT, k] is the vectore of exposure increment x_il - x_(i-1)l for patient l, expo variable k    
  # eta0 : vector all the coef of WCE
  # iWbeg, iWend : coef of the ith WCE variable is eta0[iWbeg[i]:iWend[i]]
  # ISpline_W : list of the nW integrated splines parameters for the WCE effects multiplied by eta0
  #           : thus no nead to multiply each spline coordinate by its coef
  # Spline_t0 : splines parameters for the baseline hazard multiplied by gamma0
  #           : thus no nead to multiply each spline coordinate by its coef

#  print("rateTD_gamma0alphabetaeta0")
  # spline bases for baseline hazard
    WCE <- predictSpline(Spline_t0, T, intercept=Intercept_t0, outer.ok=TRUE)
  # spline bases for each TD effect
    YT <- fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE)


    for(iW in 1:nW){
      for(iId in FirstId[iT]:iT){
        WCE <- WCE - W[iId, iW] * predictSpline(ISpline_W[[iW]], T-fromT[iId], intercept=Intercept_W[[iW]], outer.ok=TRUE)  
      }
    }
  # returned value 
  exp(WCE + YT %*% Zalphabeta[iT,])
  
}



# computes the contribution of gamma0 and eta0 (baseline & WCE)
rateTD_gamma0eta0<- function(T, iT,
                                      fromT, FirstId,
                                      gamma0, 
                                      nW, W, eta0, iWbeg, iWend,
                      Spline_t0=SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t0=TRUE,
                                      ISpline_W, Intercept_W=rep(TRUE,nW), ...){
  # compute the contribution of the time dependent variables to the rate 
  # of relative survival model for patient iT 
  # at a vector of T (useful to compute numerical integration
  # fromT : begining of the time intervals
  # toT   : end of the time intervals
  # fail  : fail = 1 if event, 0 if censored
  # FirstId : all lines in FirstId[iT]:iT of fromT, toT, fail and Zalphabeta comes from the same individual 
 # nW number of cols in W (number of WCE effects
  # W matrix of exposure INCREMENT variables W[FirstId[iT]:iT, k] is the vectore of exposure increment x_il - x_(i-1)l for patient l, expo variable k    
  # eta0 : vector all the coef of WCE
  # iWbeg, iWend : coef of the ith WCE variable is eta0[iWbeg[i]:iWend[i]]
  # ISpline_W : list of the nW integrated splines parameters for the WCE effects multiplied by eta0
  #           : thus no nead to multiply each spline coordinate by its coef
  # Spline_t0 : splines parameters for the baseline hazard multiplied by gamma0
  #           : thus no nead to multiply each spline coordinate by its coef

#  print("rateTD_gamma0eta0")
  
  # spline bases for baseline hazard
    WCE <- predictSpline(Spline_t0, T, intercept=Intercept_t0, outer.ok=TRUE)
    for(iW in 1:nW){
      for(iId in FirstId[iT]:iT){
        WCE <- WCE + W[iId, iW] * predictSpline(ISpline_W[[iW]], T-fromT[iId], intercept=Intercept_W[[iW]], outer.ok=TRUE)  
      }
    }
  # returned value 
  exp(WCE)
  
}






