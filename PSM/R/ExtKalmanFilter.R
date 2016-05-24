`ExtKalmanFilter` <-
function( phi, Model, Data, outputInternals=FALSE) {

  # Extract Data components
  Time  <- Data[["Time"]] 
  Y     <- Data[["Y"]]

  # Check for INPUT and set it
  if(is.null(Data[["U"]])) { #check if U exists.
    ModelHasInput <- FALSE
  }
  else if ( any(is.na(Data[["U"]])) ) { #check if it contains any NA
    ModelHasInput <- FALSE
  }
  else {
    ModelHasInput <- TRUE
  }
  U     <- if( !ModelHasInput) { NA } else { Data[["U"]] }

  # Set Uk
  if(ModelHasInput) {
    Uk <- U[,1,drop=FALSE]
  } else {Uk <- NA}

  DataHasDose <- "Dose" %in% names(Data)

  f  <- Model$Functions$f
  df <- Model$Functions$df
  g  <- Model$Functions$g
  dg <- Model$Functions$dg

  # Calculate Init state and initial SIG
  InitialState  <- Model$X0(Time=Time[1], phi=phi, U = Uk)
  SIG           <- Model$SIG(u=Uk,time=Time[1],phi=phi)
  Ax0 <- df(x=InitialState ,u=Uk,time=Time[1],phi=phi) 


  dimT  <- length(Time)     # Time is vector -> use length
  dimY  <- nrow(Y)          # Dimensionality of observations
  dimU  <- ifelse(ModelHasInput,nrow(U),0)
  dimX  <- nrow(InitialState)

  # P0 CTSM MathGuide page 19 (1.118) and page 8 (1.49)
  PS  <- 1.0;
  tau <- Time[2] - Time[1]

  # Dimensions Pintegral:[2*dimX 2*dimX]
  tmp <-  rbind( cbind(-Ax0 , SIG%*%t.default(SIG)) ,
                cbind( matrix(rep(0,2*dimX),nrow=dimX,ncol=dimX) , t.default(Ax0) ))*tau
  
  # Use Matrix package to compute Matrix exponential
  Pint  <- matexp(tmp)
  
  PHI0   <- t.default(Pint[(dimX+1):(2*dimX),(dimX+1):(2*dimX),drop=FALSE])
  P0    <- PS*PHI0 %*% Pint[1:dimX,(dimX+1):(2*dimX),drop=FALSE]

  #----------------------------------------------------------
  # Init matrices used in the Kalman filtering
  #----------------------------------------------------------
  #Xp <- Xf  <- array(NA,c(dimX,dimT))
  #Yp        <- array(NA,c(dimY,dimT))
  #KfGain    <- array(NA,c(dimX,dimY,dimT))
  #Pp <- Pf  <- array(NA,c(dimX,dimX,dimT))  
  #R         <- array(NA,c(dimY,dimY,dimT))

  BigNA   <- rep(NA,max(dimX,dimY)*max(dimX,dimY)*dimT)

  Xf      <- BigNA[1:(dimX*dimT)]
  dim(Xf) <- c(dimX,dimT)
  Xp <-   Xf
  Yp      <- BigNA[1:(dimY*dimT)]
  dim(Yp) <- c(dimY,dimT)
  KfGain  <- BigNA[1:(dimX*dimY*dimT)]
  dim(KfGain) <- c(dimX,dimY,dimT)
  Pf      <- BigNA[1:(dimX*dimX*dimT)]
  dim(Pf) <- c(dimX,dimX,dimT)
  Pp <- Pf
  R       <- BigNA[1:(dimY*dimY*dimT)]
  dim(R)  <- c(dimY,dimY,dimT)


  # Insert initial estimates into matrices
  Pp[,,1] <- P0
  Xp[,1]  <- InitialState
                      
  # Internal variables                          
  matXX <- matrix(0,ncol=dimX,nrow=dimX)
  DIAGDIMY  <- diag(1,dimY)

  Index <- NULL
  for (p in 1:dimX) {
    Index <- c(Index , (1:p) + dimX*(p-1) )
  }

  dSystemPred <- function(t,y,parms) {
  # Phi and Index are set in the outer enviroment
  # Only Uk needs to be updated at every observation
      Uk <- parms
      # Evaluate dX
      X  <- y[1:dimX]
      dX <- f(x=X,u=Uk,time=t,phi=phi)
      # Evaluate dP
      tmpP <- matXX
      tmpP[Index] <- y[-(1:dimX)]
      if(dimX>1) {
          tmpP[lower.tri(tmpP)] <- tmpP[upper.tri(tmpP)]
      }
      Ax <- df(x=X,u=Uk,time=t,phi=phi)
      SIGx <- Model$SIG(u=Uk,time=t,phi=phi)
      ### MathGuide (1.79)
      dP <- Ax%*%tmpP + tmpP*t.default(Ax) + SIGx%*%t.default(SIGx) 
      dP <- dP[Index]
      # Return dX and dP
      list(c( dX, dP))
    }
  
  
  ######################
  # Loop over timepoints
  ######################
  negLogLike <- 0 
  for(k in 1:dimT) {

    # Does the observation has missing values?
    ObsIndex  <- which(!is.na(Y[,k]))
    E         <- DIAGDIMY[ObsIndex,,drop=FALSE]
        
    # Set Uk
    if(ModelHasInput) {
      Uk <- U[,k,drop=FALSE]
    } else {Uk <- NA}                    
    
    
    # Y_Hat, Prediction
    Yp[ObsIndex,k] <- E %*% g(x=Xp[,k,drop=FALSE],u=Uk,time=Time[k],phi=phi)
    
    # Find the h-derivative in this point
    C <- dg(x=Xp[,k,drop=FALSE],u=Uk,time=Time[k],phi=phi)
    S <- Model$S(u=Uk,time=Time[k],phi=phi)

    # Create tmp-variable to save computation
    mattmp  <- CutThirdDim(Pp[,,k,drop=FALSE])%*%t.default(C)


    # Uncertainty on Measurement.
    R[ObsIndex,ObsIndex,k] <- E%*%C%*%mattmp%*%t.default(E) + E%*%S%*%t.default(E)


    if(length(ObsIndex)>0) { #there must be at least one obs for updating.

      ######################
      # Update
      ######################
    
      InvR    <- solve(CutThirdDim(R[ObsIndex,ObsIndex,k,drop=FALSE]))
    
      # Kalman gain
      KfGain[,ObsIndex,k] <- mattmp %*% t.default(E) %*% InvR
    
      # Updating
      e       <- Y[ObsIndex,k,drop=FALSE] - Yp[ObsIndex,k,drop=FALSE]
      KFg     <- CutThirdDim(KfGain[,ObsIndex,k,drop=FALSE])
      Xf[,k]  <- Xp[,k,drop=FALSE] + KFg%*%e
      Pf[,,k] <- Pp[,,k] - KFg %*% R[ObsIndex,ObsIndex,k] %*% t.default(KFg)

      # Add contribution to negLogLike
      tmpR       <-  CutThirdDim(R[ObsIndex,ObsIndex,k,drop=FALSE])
      logdet2piR <-  determinant.matrix(2*pi*tmpR)$modulus
      negLogLike <- negLogLike + .5*( logdet2piR + t.default(e) %*% InvR %*% e)

    } else {
      # No observations, update not available.
      Xf[,k] <- Xp[,k]
      Pf[,,k] <- Pp[,,k]
    }
    
    
    # Abort state pred if finished
    if(k==dimT) break

    
    ######################
    # Prediction
    ######################

    Xstart <- Xf[,k]
    # Add dose before starting to predict.
    if(DataHasDose) {
      # Check if dosing is occuring at this timepoint.
      if( any(Time[k]==Data$Dose$Time)) {
        idxD = which(Time[k]==Data$Dose$Time)
        # Multiple dosing a timepoint[k]
        for(cmt in 1:length(idxD)) {
          Xstart[Data$Dose$State[idxD[cmt]]] <-
            Xstart[Data$Dose$State[idxD[cmt]]] + Data$Dose$Amount[idxD[cmt]]
        }
      }
    }
   
    # Create Z combined variable
    # Upper triangle of P
    tmpP <- Pf[,,k,drop=FALSE]
    tmpP <- tmpP[Index]
    Z <- c( Xstart , tmpP)

    # Prediction of Z
    timevec <- c(Time[k], Time[k+1]) 
    ZOUT <- lsoda(y=Z, times=timevec, func=dSystemPred, parms=Uk, rtol=1e-6, atol=1e-6)
    
    # convert back to X,Pk
    Xp[,k+1] <- ZOUT[length(timevec),1+(1:dimX)] #first col is time
    
    tmpP <- matXX
    tmpP[Index] <- ZOUT[length(timevec),-(1:(dimX+1))]
    if(dimX>1) {
      tmpP[lower.tri(tmpP)] <- tmpP[upper.tri(tmpP)]
    }
    Pp[,,k+1] <- tmpP
    
  } #end loop over observations

  if(outputInternals) {
    return( list( negLogLike=negLogLike,Time=Time,Xp=Xp, Xf=Xf, Yp=Yp, KfGain=KfGain,Pf=Pf, Pp=Pp, R=R))
  } else {
    return(as.vector(negLogLike))
  }

}
