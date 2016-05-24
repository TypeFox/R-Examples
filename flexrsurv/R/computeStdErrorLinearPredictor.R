.computeStdErrorLinearPredictor_GA0B0AB<-function(GA0B0AB,
                                                  var,
                                                  Y, X0, X, Z, 
                      nT0basis,
                      Spline_t0=MSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE), Intercept_t0=TRUE,
                      ialpha0, nX0,
                      ibeta0, nX,
                      ialpha, ibeta,
                                                  nTbasis,
                      Spline_t =MSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE),
                                                  listZSplineBasis,
                      Intercept_t_NPH=rep(TRUE, nX), 
                      debug=FALSE,  ...){
  # compute std error of a linearpredictor (log rate) if the model
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
  
if ( debug) cat("  # computinfg stderr of the linear predictor: .computeStdErrorLlinearPredictor\n")

nbeta0 <- length(ibeta0)
nalpha0 <- length(ialpha0)

  # spline bases for baseline hazard
colEndTime <- ifelse(ncol(Y)==2, 1, 2)
YT0 <- evaluate(Spline_t0, Y[,colEndTime], intercept=Intercept_t0)

DesignMatrix <- cbind(YT0,
                      X0)

  # spline bases for each TD effect
    if(nX){
      # spline bases for each TD effect
      YT <- evaluate(Spline_t, Y[,colEndTime], intercept=TRUE)
      if(Intercept_t_NPH[1]){
        cbindDiagT <- YT
      }
      else {
        cbindDiagT <- YT[, -1]
      }
      
      cbindDiagX <- as.matrix(rep(1, nTbasis +  Intercept_t_NPH[1] -1 ))
      if( nX > 1 ){
        for( i in 2:nX){
          if(Intercept_t_NPH[i]){
            cbindDiagT <- cbind(cbindDiagT, YT)
          }
          else {
            cbindDiagT <- cbind(cbindDiagT, YT[, -1])
          }
          cbindDiagX <- cbindDiagX %s% as.matrix(rep(1, nTbasis + Intercept_t_NPH[i] -1))
        }
      }
      # remove first basis if
      
      DesignMatrix <- cbind(DesignMatrix,
                            (X %*%  t(cbindDiagX)) * cbindDiagT)
    }
    
if(is.null(Z)){
  nZ <- 0
  return(sqrt(diag(DesignMatrix %*% var %*% t(DesignMatrix))))
} else {
    #number of linear coef
      # spline bases for each TD effect
  YT <- evaluate(Spline_t, Y[,colEndTime], intercept=TRUE)
  nlincoef <- nT0basis + nalpha0 + nbeta0
  nZ <- Z@nZ
  alpha <- GA0B0AB[ialpha]
  beta <- GA0B0AB[ibeta]
  dGdBeta <- diag(nlincoef)
  for(i in 1:nZ){
    thenalpha <- getNBases(listZSplineBasis[[i]])-1
    cbindDiag1 <- duplicMat(diag(thenalpha), nTbasis)
    cbindDiag2 <- duplicMat(diag(nTbasis), thenalpha)
    thealpha <- alpha[Z@index[i,1]:Z@index[i,2]]
    thebeta <- beta[(i-1)*(nTbasis-1)+1:(nTbasis-1)]
    BdiagA <- duplicSumDirect(c(1, thebeta), thenalpha)
    BdiagB <- StackDiag( thealpha, nTbasis)[,-1]
    if(i == 1){
      GA <- BdiagA
      GB <- BdiagB
    }
    else {
      GA <- GA %s% BdiagA
      GA <- GB %s% BdiagB
    }
    DesignMatrix <- cbind(DesignMatrix,
                           ((Z@DM)[,Z@index[i,1]:Z@index[i,2]] %*% cbindDiag1) * (YT %*% cbindDiag2)
                           )
     dGdBeta <- dGdBeta %s% cbind(GA , GB)

    
  }

  
    DesignMatrix <- DesignMatrix %*% dGdBeta

    return(sqrt(diag(DesignMatrix %*% var %*% t(DesignMatrix))))
    
  }

}

