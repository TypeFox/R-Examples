ll_flexrsurv_fromto_beta0beta<-function(beta0beta, alpha, gamma0, alpha0,
                                        Y, X0, X, Z, 
                                        expected_rate,
                                        weights=NULL,
                                        step, Nstep, 
                               intTD=intTD_NC, intweightsfunc=intweights_CAV_SIM,
                               nT0basis,
                               Spline_t0=MSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE), Intercept_t0=TRUE, 
                               nX0,
                               ibeta0, nX, 
                               ibeta, nTbasis,
                               Spline_t =MSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE),
                               Intercept_t_NPH=FALSE,
                               debug=FALSE, ...){
  # compute log likelihood of the relative survival  model
  # rate = exp( f(t)%*%gamma + X0%*%alpha0 + X%*%beta0(t) + sum( alphai(zi)betai(t) ))
  # case where nX > 0 ie, there are NPH but LIN effects 
  # gamma : vector of coef for baseline hazard
  # alpha0 ; vector of all coefs for non time dependant variables (may contain non-loglinear terms such as spline)
  # beta0beta ; matrix of all coefs for time dependant variables  X%*%beta0(t) and alpha(Z)%*%beta(t)
  #          beta0 = beta0beta[ibeta0] and beta = beta0beta[ibeta]
                                        # beta does not contains coef for the first t-basis
  #################################################################################################################
  # alpha : vector of coef for alpha(z) for NLG and NPH
  # Y : object of class Surv with beginning and end of interval
  #
  # X0 : non-time dependante variable (may contain spline bases expended for non-loglinear terms)
  # X : log lineair but time dependante variable 
  # Z : object of class "DesignMatrixNPHNLL" time dependent variables (spline basis expended)
  # expected_rate : expected rate at event time T
  # weights : vector of weights  : LL = sum_i w_i ll_i
  # step : lag of subinterval for numerical integration fr each observation
  # Nstep : number of lag for each observation
  # intTD : function to perform numerical integration 
  # intweightfunc : function to compute weightsfor numerical integration
  # nT0basis : number of spline basis 
  #  Spline_t0, spline object for baseline hazard, with evaluate() méthod
  #  Intercept_t0=FALSE, option for evaluate, = TRUE all the basis, =FALSE all but first basis
  
  # nTbasis : number of time spline basis for NPH or NLL effects
  # nX0   : nb of PH variables dim(X0)=c(nobs, nX0)
  # nX    : nb of NPHLIN variables dim(X)=c(nobs, nX)
  #  Spline_t, spline object for time dependant effects,  with evaluate() méthod
  # Intercept_t_NPH vector of intercept option for NPH spline (=FALSE when X is NLL too, ie in case of remontet additif NLLNPH)
  #  ... not used args
  # returned value : the log liikelihood of the model
  
  if(is.null(Z)){
    nZ <- 0
  }
  else
    {
    nZ <- Z@nZ
    }

  
  # contribution of non time dependant variables
  if(nX0) PHterm <-exp(X0 %*% alpha0)
  else PHterm <- 1
  
  # contribution of time d?pendant effect
  # parenthesis are important for efficiency
  if(nZ) {
    Zalphabeta <- Z@DM %*% ( ( diag(alpha)%*% Z@signature ) %*% t(ExpandAllCoefBasis(beta0beta[ibeta], ncol=nZ,  value=1)))
    if(nX) {
      Zalphabeta <- Zalphabeta + X %*% t(ExpandCoefBasis(beta0beta[ibeta0],
                                                         ncol=nX,
                                                         splinebasis=Spline_t,
                                                         expand=!Intercept_t_NPH,
                                                         value=0))
    }
  }
  else {
    if(nX) {
      Zalphabeta <- X %*% t(ExpandCoefBasis(beta0beta[ibeta0],
                                            ncol=nX,
                                            splinebasis=Spline_t,
                                            expand=!Intercept_t_NPH,
                                            value=0))
    }
    else {
      Zalphabeta <- NULL
    }
  }
  

  if( nX+nZ ){
    NPHterm <- intTD(rateTD_gamma0alphabeta,  fromT=Y[,1], toT=Y[,2], fail=Y[,3],
                     step, Nstep, intweightsfunc=intweightsfunc, 
                     gamma0=gamma0, Zalphabeta=Zalphabeta, 
                     Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
                     Spline_t = Spline_t, Intercept_t=TRUE)
  }
  else {
    NPHterm <- intTD(rateTD_gamma0,  fromT=Y[,1], toT=Y[,2], fail=Y[,3],
                     step=step, Nstep=Nstep, intweightsfunc=intweightsfunc, 
                     gamma0=gamma0,
                     Spline_t0=Spline_t0, Intercept_t0=Intercept_t0)
  }
  # spline bases for baseline hazard at the end of intervals
  YT0 <- evaluate(Spline_t0, Y[,2], intercept=Intercept_t0)
  if(nX0){
    if(nX + nZ){
      # spline bases for each TD effect at the end of intervals
      YT <- evaluate(Spline_t, Y[,2], intercept=TRUE)
      eventterm <- ifelse(Y[,3] ,
                          log( PHterm * exp(YT0 %*% gamma0 + apply(YT * Zalphabeta, 1, sum)) + expected_rate ), 
                          0)
    } 
    else {
      eventterm <- ifelse(Y[,3] , 
                          log( PHterm * exp(YT0 %*% gamma0) + expected_rate ), 
                          0)
    }
  } 
  else {
    if(nX + nZ){
      # spline bases for each TD effect at the end of intervals
      YT <- evaluate(Spline_t, Y[,2], intercept=TRUE)
      eventterm <- ifelse(Y[,3] ,
                          log( exp(YT0 %*% gamma0 + apply(YT * Zalphabeta, 1, sum)) + expected_rate ), 
                          0)
    } else {
      eventterm <- ifelse(Y[,3] , 
                          log( exp(YT0 %*% gamma0) + expected_rate ), 
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


