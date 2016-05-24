`LinKalmanSmoother` <-
function(phi, Model, Data) {
  # Linear Continuos Kalman Filter for a state space formulation
  # Bryson-Fraizier smoother for a single subject


  # Do a forward Kalman Filter
  Obj <- LinKalmanFilter( phi , Model , Data , echo=FALSE, outputInternals=TRUE)

  

  Time  <- Data[["Time"]]
  Y     <- Data[["Y"]]
  
  # Check for INPUT and set it
  if(is.null(Data[["U"]])) { #check if U exists.
    ModelHasInput <- FALSE
  }
  else if ( sum(is.na(Data[["U"]])) >=1 ) { #check if it contains any NA
    ModelHasInput <- FALSE
  }
  else {
    ModelHasInput <- TRUE
  }
  if(ModelHasInput) {
    U <- Data[["U"]]
    Uk <- U[,1,drop=FALSE]
  } else {
    U <- NA
    Uk <- NA
  } 

  tmpM  <- Model$Matrices(phi=phi)
  matA  <- tmpM$matA
  matB  <- tmpM$matB
  matC  <- tmpM$matC
  matD  <- tmpM$matD 
  
                                        # Output prediction covariance
  S   <- Model$S(phi=phi)
  X0  <- Model$X0(Time=Time[1], phi=phi, U = Uk)
  SIG <- Model$SIG(phi=phi)


  dimN  <- length(Time)     # Time is vector -> use length
  dimY  <- nrow(Y)          # Dimensionality of observations
  dimU  <- ifelse(ModelHasInput,nrow(U),0)
  dimX  <- nrow(X0)


                                        # Init Smoothing arrays and variables
  Xs <- array(NA, c(dimX,dimN))
  Ys <- array(NA, c(dimY,dimN))
  Ps <- array(NA, c(dimX,dimX,dimN))
  lambda <- array(0.0, c(dimX))
  LAMBDA <- array(0.0, c(dimX, dimX))

  # Perform backwards filtering.
  for( i in 1:dimN) {
    
    tau <- (dimN+1)-i
    ts <- ifelse( i<dimN , Time[tau]-Time[tau-1] , Time[2]-Time[1])
    
    Adis <- matexp(matA*ts)
 
    if( any(is.na(Y[,tau,drop=FALSE]))) {  #skip update if Y_k contains any NA
      Fk <- Adis
      lambda <- t.default(Fk)%*%lambda
      LAMBDA <- t.default(Fk) %*% LAMBDA %*% Fk
      
    } else {
      KpGain <- Adis %*% CutThirdDim( Obj$KfGain[,,tau,drop=FALSE] )
      Fk <- Adis - KpGain %*% matC
      lambda <- t.default(Fk)%*%lambda +
        t.default(matC) %*% solve(Obj$R[,,tau]) %*% (Y[,tau,drop=FALSE]-Obj$Yp[,tau,drop=FALSE])
      LAMBDA <- t.default(Fk) %*% LAMBDA %*% Fk +
        t.default(matC)%*% solve(Obj$R[,,tau]) %*% matC
    }

    
    
    # Create Smooth State
    Xs[,tau] <- Obj$Xp[,tau,drop=FALSE] + CutThirdDim(Obj$Pp[,,tau,drop=FALSE]) %*% lambda

    # Create Smoothed covariance
    Ps[,,tau] <- Obj$Pp[,,tau] - Obj$Pp[,,tau] %*% LAMBDA %*% Obj$Pp[,,tau]
    Ps[,,tau] <- (Ps[,,tau] + t.default(Ps[,,tau])) / 2
    
  } # Loop over DimN

  # Create Smooth Output
  if(ModelHasInput) {
    Ys = matC%*%Xs+matD%*%U
  } else {
    Ys = matC%*%Xs
  }
  
  return( list(Time=Time, Xs=Xs, Ps=Ps, Ys = Ys, Xf=Obj$Xf, Pf=Obj$Pf, Xp=Obj$Xp, Pp=Obj$Pp,
               Yp=Obj$Yp, R=Obj$R))

}

