.computeLinearPredictor_GA0B0AB<-function(GA0B0AB, Y, X0, X, Z, 
                      nT0basis,
                      Spline_t0=MSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE), Intercept_t0=TRUE,
                      ialpha0, nX0,
                      ibeta0, nX,
                      ialpha, ibeta,                             
                      nTbasis,
                      Spline_t =MSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE),
                      Intercept_t_NPH=rep(TRUE, nX), 
                      debug=FALSE,  ...){
  # compute linearpredictor (log rate) if the model
  # rate = f(t)%*%gamma + X0%*%alpha0 + X%*%beta0(t) + sum( alphai(zi)betai(t) ))
  #################################################################################################################
  #################################################################################################################
  #  the coef of the first t-basis is constraint to 1 for nat-spline, and n-sum(other beta) if bs using expand() method
  #################################################################################################################
  #################################################################################################################
  #################################################################################################################
  # GA0B0AB ; vector of all coefs
                                        # gamma0 = GA0B0AB[1:nY0basis]
                                        # alpha0= GA0B0AB[ialpha0]
                                        # beta0= matrix(GA0B0AB[ibeta0], ncol=nX, nrow=nTbasis)
                                        # alpha= diag(GA0B0AB[ialpha])
                                        # beta= expand(matrix(GA0B0AB[ibeta], ncol=Z@nZ, nrow=nTbasis-1))
                                        # beta does not contains coef for the first t-basis
  #################################################################################################################
  # Y : object of class Surv
  # X0 : non-time dependante variable (may contain spline bases expended for non-loglinear terms)
  # X : log lineair but time dependante variable 
  # Z : object of class "DesignMatrixNPHNLL" time dependent variables (spline basis expended)
  # nT0basis : number of spline basis 
  #  Spline_t0, spline object for baseline hazard, with evaluate() méthod
  #  Intercept_t0=FALSE, option for evaluate, = TRUE all the basis, =FALSE all but first basis 
  # nTbasis : number of time spline basis for NPH or NLL effects
  # nX0   : nb of PH variables dim(X0)=c(nobs, nX0)
  # nX    : nb of NPHLIN variables dim(X)=c(nobs, nX)
  #  Spline_t, spline object for time dependant effects,  with evaluate() méthod
  # Intercept_t_NPH vector of intercept option for NPH spline (=FALSE when X is NLL too, ie in case of remontet additif NLLNPH)
  #  ... not used args
  # the function do not check the concorcance between length of parameter vectors and the number of knots and the Z.signature
  # returned value : the log liikelihood of the model
  
  if(is.null(Z)){
    nZ <- 0
  } else {
    nZ <- Z@nZ
    }

  


  # contribution of non time dependant variables
  if( nX0){
    PHterm <-X0 %*% GA0B0AB[ialpha0]
  } else {
    PHterm <- 0.0
  }
  # contribution of time d?pendant effect
  # parenthesis are important for efficiency
  if(nZ) {
    # add a row of one for the first T-basis 
    Beta <- t(ExpandAllCoefBasis(GA0B0AB[ibeta], ncol=nZ,  value=1))
    # parenthesis important for speed ?
    Zalphabeta <- Z@DM %*%( diag(GA0B0AB[ialpha]) %*% Z@signature  %*% Beta )
    if(nX) {
    # add a row of 0 for the first T-basis when !Intercept_T_NPH
      Zalphabeta <- Zalphabeta + X %*% t(ExpandCoefBasis(GA0B0AB[ibeta0],
                                                         ncol=nX,
                                                         splinebasis=Spline_t,
                                                         expand=!Intercept_t_NPH,
                                                         value=0))
    }
  } else {
    if(nX) {
      Zalphabeta <- X %*% t(ExpandCoefBasis(GA0B0AB[ibeta0],
                                            ncol=nX,
                                            splinebasis=Spline_t,
                                            expand=!Intercept_t_NPH,
                                            value=0))
    }
  }
  
  # spline bases for baseline hazard
  colEndTime <- ifelse(ncol(Y)==2, 1, 2)
    YT0 <- evaluate(Spline_t0, Y[,colEndTime], intercept=Intercept_t0)
    
  # spline bases for each TD effect
    if(nX + nZ){
      # spline bases for each TD effect
      YT <- evaluate(Spline_t, Y[,colEndTime], intercept=TRUE)
      linpred <- YT0 %*% GA0B0AB[1:nT0basis] + PHterm + apply(YT * Zalphabeta, 1, sum)
    } else {
      linpred <- YT0 %*% GA0B0AB[1:nT0basis] + PHterm 
    }


linpred
}






