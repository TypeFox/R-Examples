gr_ll_flexrsurv_fromto_GA0B0AB<-function(GA0B0AB, Y, X0, X, Z, 
                                         expected_rate,
                                         weights=NULL,
                                         step, Nstep, 
                                         intTD=intTD_NC, intweightsfunc=intweights_CAV_SIM,
                                         intTD_base=intTD_base_NC,
                      nT0basis,
                      Spline_t0=MSplineBasis(knots=NULL, degree=3,   keep.duplicates=TRUE), Intercept_t0=TRUE,
                      ialpha0, nX0,
                      ibeta0, nX,
                      ialpha, ibeta,                             
                      nTbasis,
                      Spline_t =MSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE),
                      Intercept_t_NPH=rep(TRUE, nX),
                                                     debug.gr=FALSE,  ...){
  # compute gradient of the log likelihood of the relatice survival model
  # rate = exp( f(t)%*%gamma + X0%*%alpha0 + X%*%beta0(t) + sum( alphai(zi)betai(t) ))
  #################################################################################################################
  #################################################################################################################
  #  the coef of the first t-basis is constraint to 1 for nat-spline, and n-sum(other beta) if BS using expand() method
  #################################################################################################################
  #################################################################################################################
  #################################################################################################################
  # GA0B0AB ; vector of all coefs
                                        # gamma0 = GA0B0AB[1:nY0basis]
                                        # alpha0= GA0B0AB[ialpha0]
                                        # beta0= matrix(GA0B0AB[ibeta0], ncol=nX, nrow=nTbasis)
                                        # alpha= diag(GA0B0AB[ialpha])
                                        # beta= expand(matrix(GA0B0AB[ibeta], ncol=nZ, nrow=nTbasis-1))
                                        # beta does not contains coef for the first t-basis
  #################################################################################################################
  # Y : object of class Surv with beginning and end of interval
  #
  # X0 : non-time dependante variable (may contain spline bases expended for non-loglinear terms)
  # X : log lineair but time dependante variable 
  # Z : objesct of class DeSignMatrixLPHNLL of time dépendent variables (spline basis expended)
  # expected_rate : expected rate at event time T
  # weights : weights  : LL = sum_i w_i ll_i
  # step : lag of subinterval for numerical integration fr each observation
  # Nstep : number of lag for each observation
  # intTD : function to perform numerical integration 
  # intweightfunc : function to compute weightsfor numerical integration
  # Knots_t0=NULL,Intercept_t0=FALSE, degree_t0=3, Boundary.knots_t0 time spline parameters for baseline hazard
  # Knots_t=NULL,Intercept_t=FALSE, degree_t0=, Boundary.knots_t  time spline parameters for time-dependant effects (same basis for each TD variable)
  # nT0basis : number of spline basis for NPHLIN effects
  # nX0   : nb of PH variables dim(X0)=c(nobs, nX0)
  # nX    : nb of NPHLIN variables dim(X)=c(nobs, nX)
  # nTbasis : number of time spline basis
  #  ... not used args
  # the function do not check the concorcance between length of parameter vectors and the number of knots and the Z@signature
  # returned value : the log liikelihood of the model
  
  if (debug.gr) cat("# computing gradient: gr_ll_flexrsurv_GA0B0AB\n")

  if(is.null(Z)){
    nZ <- 0
  }
  else {
    nZ <- Z@nZ
  }
 



  # contribution of non time dependant variables
  if( nX0){
    PHterm <-exp(X0 %*% GA0B0AB[ialpha0])
  }
  else PHterm <- 1
  # contribution of time d?pendant effect
  # parenthesis are important for efficiency
  if(nZ) {
    # add a row for the first basis
    tBeta <- t(ExpandAllCoefBasis(GA0B0AB[ibeta], ncol=nZ,  value=1))
    # Zalpha est la matrice des alpha(Z)
    # parenthesis important for speed ?
    Zalpha <- Z@DM %*%( diag(GA0B0AB[ialpha]) %*% Z@signature )
    Zalphabeta <- Zalpha  %*% tBeta 
    if(nX) {
      Zalphabeta <- Zalphabeta + X %*% t(ExpandCoefBasis(GA0B0AB[ibeta0],
                                                         ncol=nX,
                                                         splinebasis=Spline_t,
                                                         expand=!Intercept_t_NPH,
                                                         value=0))
    }
  }
  else {
    if(nX) {
      Zalphabeta <- X %*% t(ExpandCoefBasis(GA0B0AB[ibeta0],
                                            ncol=nX,
                                            splinebasis=Spline_t,
                                            expand=!Intercept_t_NPH,
                                            value=0))
    }
  }
  
  if(nX + nZ) {
    NPHterm <- intTD(rateTD_gamma0alphabeta, fromT=Y[,1], toT=Y[,2], fail=Y[,3],
                     step=step, Nstep=Nstep, intweightsfunc=intweightsfunc, 
                     gamma0=GA0B0AB[1:nT0basis], Zalphabeta=Zalphabeta, 
                     Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
                     Spline_t = Spline_t, Intercept_t=TRUE)
    Intb0 <-  intTD_base(func=rateTD_gamma0alphabeta, fromT=Y[,1], toT=Y[,2], fail=Y[,3],
                         Spline=Spline_t0,
                         step=step, Nstep=Nstep, intweightsfunc=intweightsfunc, 
                         gamma0=GA0B0AB[1:nT0basis], Zalphabeta=Zalphabeta, 
                         Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
                         Spline_t = Spline_t, Intercept_t=TRUE,
                         debug=debug.gr)
    if( identical(Spline_t0, Spline_t)){
      Intb <- Intb0
    }
    else {
      Intb <-  intTD_base(func=rateTD_gamma0alphabeta, fromT=Y[,1], toT=Y[,2], fail=Y[,3],
                          Spline=Spline_t,
                          step=step, Nstep=Nstep, intweightsfunc=intweightsfunc,
                          gamma0=GA0B0AB[1:nT0basis], Zalphabeta=Zalphabeta, 
                          Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
                          Spline_t = Spline_t, Intercept_t=TRUE)
    }
    if(!Intercept_t0){
     Intb0<- Intb0[,-1]
    }
    indx_without_intercept <- 2:getNBases(Spline_t)
    
    YT0 <- fevaluate(Spline_t0, Y[,2], intercept=Intercept_t0)
    YT <- fevaluate(Spline_t, Y[,2], intercept=TRUE)
    RatePred <- ifelse(Y[,3] ,
                       PHterm * exp(YT0 %*% GA0B0AB[1:nT0basis] + apply(YT * Zalphabeta, 1, sum)),
                       0)
  }  # weights : weights  : LL = sum_i w_i ll_i

  else {
    NPHterm <- intTD(rateTD_gamma0, fromT=Y[,1], toT=Y[,2], fail=Y[,3],
                     step=step, Nstep=Nstep, intweightsfunc=intweightsfunc, 
                     gamma0=GA0B0AB[1:nT0basis],
                     Spline_t0=Spline_t0, Intercept_t0=Intercept_t0)
    Intb0 <-  intTD_base(func=rateTD_gamma0, fromT=Y[,1], toT=Y[,2], fail=Y[,3],
                         Spline=Spline_t0,
                         step=step, Nstep=Nstep, intweightsfunc=intweightsfunc, 
                         gamma0=GA0B0AB[1:nT0basis], 
                         Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
                         debug=debug.gr)
    if(!Intercept_t0){
     Intb0<- Intb0[,-1]
    }
    Intb <- NULL
    YT0 <- fevaluate(Spline_t0, Y[,2], intercept=Intercept_t0)
    YT <- NULL
    RatePred <- ifelse(Y[,3] ,
                       PHterm * exp(YT0 %*% GA0B0AB[1:nT0basis]) ,
                       0)
                         
  }

  F <- ifelse(Y[,3] ,
              RatePred/(RatePred + expected_rate ), 
              0)
  
  
  if(nX + nZ) {
    if(nX0>0) {
      Intb <- Intb * c(PHterm)
  }
    IntbF <- YT*F - Intb
  }
  else {
    IntbF <- NULL
  }
  Intb0 <- Intb0 * c(PHterm)
  

#####################################################################"
# now computes the mean score


  if (!is.null(weights)) {
# dldgamma0
    dLdgamma0 <- crossprod( YT0 * F - Intb0 , weights)

    if (nX0) {
      dLdalpha0 <- crossprod(X0 , (F - PHterm * NPHterm) * weights )
    }
    else {
      dLdalpha0 <- NULL
    }

    if (nX){
#  traiter les Intercept_t_NPH
      dLdbeta0 <- NULL
      for(i in 1:nX){
        if ( Intercept_t_NPH[i] ){
          dLdbeta0 <- c(dLdbeta0,  crossprod(X[,i] ,  IntbF * weights))
        }
        else {
          dLdbeta0 <- c(dLdbeta0, crossprod(X[,i] ,  IntbF[,indx_without_intercept] * weights))
        }
      }
    }
    else {
      dLdbeta0 <- NULL
    }
    
    if (nZ) { 
      baseIntbF <- IntbF  %*% t(tBeta)
      dLdalpha <- rep(0,getNparam(Z) )
      indZ <- getIndex(Z)
      
      for(iZ in 1:nZ){
        if ( debug.gr > 200 ){
          
        }
        dLdalpha[indZ[iZ,1]:indZ[iZ,2]] <- crossprod(Z@DM[,indZ[iZ,1]:indZ[iZ,2]], baseIntbF[,iZ] * weights )
      }
      dLdbeta <- c(crossprod((IntbF[,-1, drop=FALSE]),Zalpha * weights))
    }
    else {
      dLdalpha <- NULL
      dLdbeta <- NULL
    }
  } # end weights!=NULL
  else {
# d<dgamma0
    dLdgamma0 <- apply(   YT0 * F - Intb0 , 2, sum)

    if (nX0) {
      dLdalpha0 <- crossprod(X0 , F - PHterm* NPHterm )
    }
    else {
      dLdalpha0 <- NULL
    }

    if (nX){
#  traiter les Intercept_t_NPH
      dLdbeta0 <- NULL
      for(i in 1:nX){
        if ( Intercept_t_NPH[i] ){
          dLdbeta0 <- c(dLdbeta0,  crossprod(X[,i] ,  IntbF))
        }
        else {
          dLdbeta0 <- c(dLdbeta0, crossprod(X[,i] ,  IntbF[,indx_without_intercept]))
        }
      }
    }
    else {
      dLdbeta0 <- NULL
    }
    
    if (nZ) { 
      baseIntbF <- IntbF  %*% t(tBeta)
      dLdalpha <- rep(0,getNparam(Z) )
      indZ <- getIndex(Z)
      
      for(iZ in 1:nZ){
        dLdalpha[indZ[iZ,1]:indZ[iZ,2]] <- crossprod(Z@DM[,indZ[iZ,1]:indZ[iZ,2]], baseIntbF[,iZ] )
      }
      dLdbeta <- c(crossprod((IntbF[,-1, drop=FALSE]),Zalpha ))
    }
    else {
      dLdalpha <- NULL
      dLdbeta <- NULL
    }
  } # end weights==NULL

  


  rep <- c(dLdgamma0,          
           dLdalpha0,          
           dLdbeta0,          
           dLdalpha,          
           dLdbeta )

  if(debug.gr){
    attr(rep, "intb0") <- Intb0
    attr(rep, "intb") <- Intb
    attr(rep, "intbF") <- IntbF
    attr(rep, "F") <- F 
    attr(rep, "YT0") <- YT0 
    attr(rep, "YT") <- YT
    attr(rep, "RatePred") <-  RatePred
    if(debug.gr > 1000){
      cat("grad value and parameters :", "\n") 
      print(cbind(    c(dLdgamma0,          
                        dLdalpha0,          
                        dLdbeta0,          
                        dLdalpha,          
                        dLdbeta ), GA0B0AB))
      
    }
  }
  
  rep
  
}
