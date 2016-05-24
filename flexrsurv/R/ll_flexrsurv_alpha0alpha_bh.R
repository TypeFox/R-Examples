ll_flexrsurv_alpha0alpha_bh<-function(alpha0alpha, beta0, beta, gamma0,
                                   Y, X0, X, Z, 
                                   expected_rate,
                                   weights=NULL,
                                   step, Nstep, 
                                 intTD=intTD_NC, intweightsfunc=intweights_CAV_SIM,
                                 nT0basis,
                                 Spline_t0=MSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE), Intercept_t0=TRUE,
                                 ialpha0, nX0,
                                 nX, 
                                 ialpha,
                                 nTbasis,
                                 Spline_t =MSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE),
                                 Intercept_t_NPH=rep(TRUE, nX), 
                                 debug=FALSE,  ...){
  # compute log likelihood of the relative survival model
  # rate = exp( f(t)%*%gamma + X0%*%alpha0 + X%*%beta0(t) + sum( alphai(zi)betai(t) ))
  # case where nX0 > 0 ie, there are PH effects 
  # gamma : vector of coef for baseline hazard
  # alpha0alpa ; vector of all coefs for non time dependant variables (may contain non-loglinear terms 
  #              such as spline) AND alpah(Z) for NLG and NPH
                                        # alpha0= alpha0alpha[ialpha0]
                                        # alpha= alpha0alpah[ialpha]
  # beta0 ; matrix of all coefs for log-linear but  time dependant variables  X%*%beta0(t)
  # beta  : matrix of coefs for beta(t) nTbasis * nTDvars for NLG and NPH
  # Y : object of class Surv
  # X0 : non-time dependante variable (may contain spline bases expended for non-loglinear terms)
  # X : log lineair but time dependante variable 
  # Z : object of class "DesignMatrixNPHNLL" time dependent variables (spline basis expended)
  # expected_rate : expected rate at event time T
  # weights : vector of weights  : LL = sum_i w_i ll_i
  # step : lag of subinterval for numerical integration fr each observation
  # Nstep : number of lag for each observation
  # intTD : function to perform numerical integration 
  # intweightfunc : function to compute weightsfor numerical integration
  #  Spline_t0, spline object for baseline hazard, with evaluate() méthod
  #  Intercept_t0=FALSE, option for evaluate, = TRUE all the basis, =FALSE all but first basis 
  # nT0basis : number of spline basis for NPHLIN effects
  # nX0   : nb of PH variables dim(X0)=c(nobs, nX0)
  # nX    : nb of NPHLIN variables dim(X)=c(nobs, nX)
  # nTbasis : number of time spline basis
  #  Spline_t, spline object for time dependant effects,  with evaluate() méthod
  # Intercept_t_NPH vector of intercept option for NPH spline (=FALSE when X is NLL too, ie in case of remontet additif NLLNPH)
  #  ... not used args
  # the function do not check the concorcance between length of parameter vectors and the number of knots and the Z.signature
  # returned value : the log liikelihood of the model

  if(is.null(Z)){
    nZ <- 0
  }
  else
    {
    nZ <- Z@nZ
    }
  

  # contribution of non time dependant variables
  if( nX0){
    PHterm <-exp(X0 %*% alpha0alpha[ialpha0])
  }
  else PHterm <- 1
  # contribution of time d?pendant effect
  # parenthesis are important for efficiency
  if(nZ) {
    Zalphabeta <- Z@DM %*%( diag(alpha0alpha[ialpha]) %*% Z@signature  %*% t(ExpandAllCoefBasis(beta, ncol=nZ,  value=1)) )
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
    NPHterm <- intTD(rateTD_bh_alphabeta, Y[,1], fail=Y[,2],
                     step, Nstep, intweightsfunc=intweightsfunc, 
                     gamma0=gamma0, Zalphabeta=Zalphabeta, 
                     Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
                     Spline_t = Spline_t, Intercept_t=TRUE)
  }
  else {
    NPHterm <- intTD(rateTD_bh, Y[,1], fail=Y[,2],
                     step=step, Nstep=Nstep, intweightsfunc=intweightsfunc, 
                     gamma0=gamma0,
                     Spline_t0=Spline_t0, Intercept_t0=Intercept_t0)
  }
  # spline bases for baseline hazard
  YT0 <- evaluate(Spline_t0, Y[,1], intercept=Intercept_t0)
  # spline bases for each TD effect
    if(nX + nZ){
      # spline bases for each TD effect
      YT <- evaluate(Spline_t, Y[,1], intercept=TRUE)
      eventterm <- ifelse(Y[,2] ,
                          log( PHterm * (YT0 %*% gamma0) * exp(apply(YT * Zalphabeta, 1, sum)) + expected_rate ), 
                          0)
    } 
    else {
      eventterm <- ifelse(Y[,2] , 
                          log( PHterm * (YT0 %*% gamma0) + expected_rate ), 
                          0)
    }

  if (!is.null(weights)) {
    ret <- crossprod(eventterm - PHterm * NPHterm , weights)
  }
  else {
    ret <- sum( eventterm - PHterm * NPHterm )
  }


ret
}






