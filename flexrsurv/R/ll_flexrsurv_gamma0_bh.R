ll_flexrsurv_gamma0_bh<-function(gamma0, alpha0, beta0, alpha, beta, Y, X0, X, Z, 
                                 expected_rate,
                                 weights=NULL,
                                 step, Nstep,
                            intTD=intTD_NC, intweightsfunc=intweights_CAV_SIM,
                            Spline_t0=MSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE),
                              Intercept_t0=FALSE,
                            nT0basis,
                            nX0,
                            nX, 
                            nTbasis,
                            Spline_t =MSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE),
                            Intercept_t_NPH=rep(TRUE, nX),
                            debug=FALSE, ...){
  # compute log likelihood of the relative survival model model
  # rate = (f(t)%*%gamma)  * exp(  X0%*%alpha0 + X%*%beta0(t) + sum( alphai(zi)betai(t) ))
  # gamma0 : vector of coef for baseline hazard
  # alpha0 ; vector of all coefs for non time dependant variables (may contain non-loglinear terms such as spline)
  # beta0 ; matrix of all coefs for log-linear but  time dependant variables  X%*%beta0(t)
  # beta  : matrix of coefs for beta(t) nTbasis * nTDvars for NLG and NPH
  # alpha : vector of coef for alpha(z) for NLG and NPH
  # Y : object of class Surv
  # X0 : non-time dependante variable (may contain spline bases expended for non-loglinear terms)
  # X : log lineair but time dependante variable 
  # Z : object of class "DesignMatrixNPHNLL" time dependent variables (spline basis expended)
  # expected_rate : expected rate at event time T
  # weights : vector of weights  : LL = sum_i w_i ll_i
  # Step : object of class "NCLagParam" or "GLMLagParam"
  # intTD : function to perform numerical integration 
  # intweightfunc : function to compute weightsfor numerical integration
  #  Spline_t0, spline object for baseline hazard, with evaluate() méthod
  #  Intercept_t0=FALSE, option for evaluate, = TRUE all the basis, =FALSE all but first basis 
  #  Spline_t, spline object for time dependant effects,  with evaluate() méthod
  # Intercept_t_NPH vector of intercept option for NPH spline (=FALSE when X is NLL too, ie in case of remontet additif NLLNPH)
  # nT0basis : number of spline basis for NPHLIN effects
  # nX0   : nb of PH variables dim(X0)=c(nobs, nX0)
  # nX    : nb of NPHLIN variables dim(X)=c(nobs, nX)
  # nTbasis : number of time spline basis
  #  ... not used args
  # the function do not check the concorcance between length of parameter vectors and the number of knots and the Z.signature
  # returned value : the log liikelihood of the model
  
  if(is.null(Z)){
    nZ <- 0
  }
  else {
    nZ <- Z@nZ
    }


  # contribution of non time dependant variables
  if(nX0) PHterm <-exp(X0 %*% alpha0)
  else PHterm <- 1
  # contribution of time d?pendant effect
  # parenthesis are important for efficiency
  if(nZ) {
    Zalphabeta <- Z@DM %*%( diag(alpha) %*% Z@signature  %*% t(ExpandAllCoefBasis(beta, ncol=nZ,  value=1)) )
    if(nX) {
      Zalphabeta <- Zalphabeta +  X %*% t(ExpandCoefBasis(beta0,
                                                          ncol=nX,
                                                          splinebasis=Spline_t,
                                                          expand=!Intercept_t_NPH,
                                                          value=0))
    }
  }
  else {
    if(nX) {
      Zalphabeta <-  X %*% t(ExpandCoefBasis(beta0,
                                             ncol=nX,
                                             splinebasis=Spline_t,
                                             expand=!Intercept_t_NPH,
                                             value=0))
    }
    else {
      Zalphabeta <- NULL
    }
  }
  
  if(nX + nZ) {
    NPHterm <- intTD(func=rateTD_bh_alphabeta, T=Y[,1], fail=Y[,2],
                     step=step, Nstep=Nstep,  intweightsfunc=intweightsfunc, 
                     gamma0=gamma0, Zalphabeta=Zalphabeta, 
                     Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
                     Spline_t = Spline_t, Intercept_t=TRUE,
                     debug=debug)
  }
  else {
    NPHterm <- intTD(func=rateTD_bh, Y[,1], fail=Y[,2],
                     step=step, Nstep=Nstep, intweightsfunc=intweightsfunc, 
                     gamma0=gamma0,
                     Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
                     debug=debug)
  }
# spline bases for baseline hazard
  YT0 <- evaluate(Spline_t0, Y[,1], intercept=Intercept_t0)
  # spline bases for each TD effect
  if(nX0){
    if(nX + nZ){
      # spline bases for each TD effect
      YT <- evaluate(Spline_t, Y[,1], intercept=TRUE)
       if ( debug > 1000){
       }
      eventterm <- ifelse(Y[,2] ,
                          log( PHterm * (YT0 %*% gamma0) * exp(apply(YT * Zalphabeta, 1, sum)) + expected_rate ), 
                          0)
    } 
    else { 
      eventterm <- ifelse(Y[,2] , 
                          log( PHterm * (YT0 %*% gamma0) + expected_rate ), 
                          0)
    }
  } 
  else {
    if(nX + nZ){
      # spline bases for each TD effect
      YT <- evaluate(Spline_t, Y[,1], intercept=TRUE)
       if ( debug > 1000){
       }
      eventterm <- ifelse(Y[,2] ,
                          log( (YT0 %*% gamma0) * exp(apply(YT * Zalphabeta, 1, sum)) + expected_rate ), 
                          0)
    } else {
      eventterm <- ifelse(Y[,2] , 
                          log( (YT0 %*% gamma0) + expected_rate ), 
                          0)
    }
  }


  if (!is.null(weights)) {
    if( nX0){
      ret <- crossprod(eventterm - PHterm * NPHterm , weights)
    } else {
      ret <- crossprod(eventterm - NPHterm , weights)
    }
  }
  else {
    if( nX0){
      ret <- sum( eventterm - PHterm * NPHterm )
    } else {
      ret <- sum( eventterm - NPHterm )
    }
  }

ret 
}
