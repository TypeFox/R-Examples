`ExtKalmanSmoother` <-
function(phi, Model, Data) {
  # Smoothing based on CTSM MathGuide www.imm.dtu.dk/ctsm


  # Do forward EKF
  o <- ExtKalmanFilter(phi=phi, Model=Model, Data=Data, outputInternals=TRUE)

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

  DataHasDose <- "Dose" %in% names(Data)

  f  <- Model$Functions$f
  df <- Model$Functions$df
  g  <- Model$Functions$g
  dg <- Model$Functions$dg

  dimT <- length(Time)     # Time is vector -> use length
  dimX <- dim(o$Xf)[1]
  dimY  <- nrow(Y)          # Dimensionality of observations

  #Set Variables
  Ys <- array(0,c(dimY,dimT))
  sf <- array(0,c(dimX,dimT))
  sPred <- array(0,c(dimX,dimT))
  PbPred <- array(0,c(dimX,dimX,dimT))
  PbFil <- array(0,c(dimX,dimX,dimT))

  tmpXb <- sPred
  
  matXX <- matrix(NA,ncol=dimX,nrow=dimX)
  
  #Initial Conditions
  #Both S and Pback starts as Zero pr. defenition

  #Calculate indexes used to extract upperpart of covariance matrix.

  Index <- NULL
  for (p in 1:dimX) {
    Index <- c(Index , (1:p) + dimX*(p-1) )
  }


  # Prediction of Z
  # The function dSmoothPred implements the equations
  # (1.176) and (1.178) in the CTSM MathGuide
  dSmoothPred <- function(t,y,parms) {
    #Calculate foreward state from subsampling
    X1<- o$Xp[,tau-1,drop=FALSE]
    X2 <- o$Xp[,tau,drop=FALSE]

    if(t==Time[tau]) {
      XpI <- X2
    } else if (t==Time[tau-1]) {
      XpI <- X1
    } else {
      XpI <- X1 + (X2-X1)/(Time[tau]-Time[tau-1])*(t-Time[tau-1])
    }
    
    s <- y[1:dimX]

    #only Upper part given - convert into a real covariance matrix
    Pback <- matXX
    Pback[Index] <- y[-(1:dimX)]
    if(dimX>1) {
      Pback[lower.tri(Pback)] <- Pback[upper.tri(Pback)]
    }
      
    # Calculate variables
    Atau <- df(x=XpI,u=parms,time=t,phi=phi)
    SIG <- Model$SIG(u=parms,time=t,phi=phi)
    fx <- f(x=XpI,u=parms,time=t,phi=phi)
    
    # s Prediction (1.176 + 177)
    ds <- t(Atau)%*%s-Pback%*%SIG%*%t(SIG)%*%s - Pback%*%(fx-Atau%*%XpI)

    #------------------------------------------------
    #% Pback Prediction 
    #% (1.178)
    dP <- Pback%*%Atau + t(Atau)%*%Pback - Pback%*%SIG%*%t(SIG)%*%Pback
    
    #%Since Pk i symmetric it is only neccesary to Return Upper Part
    dP <- dP[Index]

    #%Collect results and return result
    list(-1*c(ds, dP))
  }

  
  # Loop over the length of the dataset
  for (k in 1:dimT) {
    
    # Tau runs from N -> 1
    tau <- (dimT+1) - k

    # Set Uk
    if(ModelHasInput) {
      Uk <- U[,k,drop=FALSE]
    } else {Uk <- NA}                    
    
    # Calculate variables
    # Formal (1.179)
    Ck <- dg(x=o$Xp[,tau,drop=FALSE],u=Uk,time=Time[tau],phi=phi)
    S <- Model$S(u=Uk,time=Time[tau],phi=phi)

    if(!any(is.na(Y[,tau,drop=FALSE]))) { #skip update if Y_k contains any NA
      # Innovations
      e <- Y[,tau,drop=FALSE] - g(x=o$Xp[,tau,drop=FALSE],u=Uk,time=Time[tau],phi=phi)

      # Updates
      # Formel (1.174) og (1.175)
      sf[,tau] <- sPred[,tau,drop=FALSE] +
        t(Ck)%*%solve(S)%*%( e + Ck%*%o$Xp[,tau,drop=FALSE] )
      PbFil[,,tau] <- PbPred[,,tau] + t(Ck)%*%solve(S)%*%Ck
    } else {
      sf[,tau] <- sPred[,tau,drop=FALSE]
      PbFil[,,tau] <- PbPred[,,tau]
    }
    
    
    #Do not predict for the last observation
    if(k==dimT) break

    
    #Prediction equations
    #Create Z combined variable

    # Extract Upper part of covariance matrix
    tmpP = PbFil[,,tau]
    tmpP = tmpP[Index]

    # Combine with s
    Z = c(sf[,tau] , tmpP)

    
    #function dZ = dSmoothPred(t,Z,Uk,theta,fktList,dimX,Index,Xf);
    #    ZOUT <- ode15s(@dSmoothPred,[T(tau) T(tau-1)], ...  
    #       Z ,[],Uk,theta,fktList,dimX,Index,o.sub_foreward(tau-1));
    timevec <- c(Time[tau],Time[tau-1]) 
    ZOUT <- lsoda(y=Z, times=timevec, func=dSmoothPred, parms=Uk, rtol=1e-6, atol=1e-6)
    
    # convert back to state and PbPred covariance matrix
    # %Store as predictions to timepoint tau-1
    
    #Store s
    sPred[,tau-1] <- ZOUT[length(timevec),(1:dimX)+1]

 
    #Store PbPred
    tmpP <- matXX
    tmpP[Index] <- ZOUT[length(timevec),-(1:(dimX+1))]
    if(dimX>1) {
      tmpP[lower.tri(tmpP)] <- tmpP[upper.tri(tmpP)]
    }
    PbPred[,,tau-1] <- tmpP


    
 
    if(DataHasDose ) {
      # Check if dosing is occuring at this timepoint.
      if( any(Time[tau-1]==Data$Dose$Time)) {
        idxD = which(Time[tau-1]==Data$Dose$Time)
        # Multiple dosing a timepoint[k]
        for(cmt in 1:length(idxD))
          sPred[,tau-1] <- sPred[,tau-1] -
            tmpP[,Data$Dose$State[idxD[cmt]]]*Data$Dose$Amount[idxD[cmt]]
      }
    }
    

  }


  #Create the Smooth output
  XSmooth <- array(0,c(dimX,dimT))
  PSmooth <- array(0,c(dimX,dimX,dimT))

  for (k in 1:dimT) {
    InvPf <- solve(o$Pf[,,k])
    
    #Mathguide (1.181)
    PSmooth[,,k] <- solve( InvPf + PbPred[,,k] )
    
    #Mathguide (1.180)
    XSmooth[,k] <- PSmooth[,,k]%*%( InvPf %*% o$Xf[,k] + sPred[,k,drop=FALSE])

    # Create Smooth Output
    if(ModelHasInput) {
      Ys[,k] <- g(x=XSmooth[,k,drop=FALSE],u=U[,k,drop=FALSE],time=Time[k],phi=phi)
    } else {
      Ys[,k] <- g(x=XSmooth[,k,drop=FALSE],u=NA,time=Time[k],phi=phi)
    }

  }

  
  
  return( list(Time=Time, Xs=XSmooth, Ps=PSmooth, Ys = Ys, Xf=o$Xf,
               Pf=o$Pf, Xp=o$Xp, Pp=o$Pp,
               Yp=o$Yp, R=o$R))

}

