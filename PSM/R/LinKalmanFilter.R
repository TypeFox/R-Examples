`LinKalmanFilter` <-
function( phi , Model , Data , echo=FALSE, outputInternals=FALSE,fast=TRUE) {
# Linear Continous Kalman Filter for a State Space formulation
#     Input: 
#     Model           type: list
#       $Matrices     function  input:  phi,u
#       $X0           function  input:  phi,u,t
#       $SIG          function  input:  phi,u,t
#
#     Data            type: list
#       $Time         matrix [1 dimT]
#       $Y            matrix [dimY dimT]
#       $U            matrix [dimU dimT]
#
#     phi             Parameter vector
#     echo            TRUE or FALSE   Display informaiton during execution
#
#     Problem definition
#     dx = (Ax + Bu)dt + SIG(phi,u,t) db
#     y = Cx + Du + e             e ~ N(0,S(phi,u,t))
#
#     Developer Notes:
#
#     Important dimensions
#       dimT    number of timepoints (observations)
#       dimX    number of states
#       dimU    dimension of input
#       dimY    output dimension
#
#     A:[dimX dimX] ; B:[dimX dimU] ; SIG:[dimX dimX]
#     C:[dimY dimX] ; D:[dimY dimU] ; S:[dimY,dimY]
#


  echo = FALSE

    # Print current value of phi .. Handy in minimizations.
  if(echo) cat("phi:" , paste(round(as.double(phi),2) , "\t"))

    # Extract Data components
    Time  <- Data[["Time"]] 
    Y     <- Data[["Y"]]

    # Check for INPUT and set it
  if(is.null(Data[["U"]])) { #check if U exists.
    ModelHasInput <- FALSE
  }  else if ( any(is.na(Data[["U"]])) ) { #check if it contains any NA
    ModelHasInput <- FALSE
  }  else {
    ModelHasInput <- TRUE
  }
  U     <- if( !ModelHasInput) { NA } else { Data[["U"]] }

    # Set Uk
    if(ModelHasInput) {
      Uk <- U[,1,drop=FALSE]
    } else {Uk <- NA}

    DataHasDose <- "Dose" %in% names(Data)
     
    tmp   <- Model$Matrices(phi=phi)
    matA  <- tmp$matA
    matB  <- tmp$matB
    matC  <- tmp$matC
    matD  <- tmp$matD 
  
  # Calculate Init state and initial SIG
  InitialState  <- Model$X0(Time=Time[1], phi=phi, U = Uk)
  SIG           <- Model$SIG(phi=phi)
  S             <- Model$S(phi=phi)
  
  dimT  <- length(Time)     # Time is vector -> use length
  dimY  <- nrow(Y)          # Dimensionality of observations
  dimU  <- ifelse(ModelHasInput,nrow(U),0)
  dimX  <- nrow(InitialState)

    # Is A singular ???
    rankA <- qr(matA)$rank
    singA <- (rankA<dimX)
    
  if(echo) {
      cat("\n *** \ndimX= " ,dimX, "\n")
      cat("dimY= " ,dimY, "\n")  
      cat("dimU= " ,dimU, "\n")
      cat("dimT= " ,dimT, "\n")
      cat("rankA= " , rankA , "\n")
      cat("singA= " , singA , "\n")
      }
      


  # Use compiled Fortran Code
  if(fast && !singA) {
  
    if(echo) cat("Preparing compiled code \n")

    #Currently only nonsingular A
    if(singA) {
      stop("Singular A not currently available in compiled code")
    }

    # Missing observations change NA -> BIG M
    Y[is.na(Y)] <- 1E300

    # Pseudo input - Fortran code assumes input
    if(!ModelHasInput) {
      if(echo) {cat("CREATING PSEUDO INPUT U \n")}
      dimU <- 1
      U <- array( 0 , c(dimU,dimT))
      matB <- array(0 , c(dimX,dimU))
      matD <- array(0 , c(dimY,dimU))      
    }
    
    # DOSING
    if(DataHasDose) {
      if(echo) {cat("USING ORIGINAL DOSING SCHEME \n")}
      DoseN     <- length(Data$Dose$Time)
      DoseTime  <- Data$Dose$Time
      DoseState <- Data$Dose$State
      DoseAmt   <- Data$Dose$Amount
    } else {
      # Create Pseudo Dose with Amount 0
      if(echo) {cat("CREATING PSEUDO DOSE \n")}
      DoseN     <- 1
      DoseTime  <- Time[1]
      DoseState <- 1
      DoseAmt   <- 0
    }
        
    # Fortran Internals - Initialize
    if(echo) {cat("CREATING FORTRAN INTERNALS \n")}
    LL <-  -1.0
    INFO <- -1
    YP <- rep(-1 , dimY*dimT)
    XF <- rep(-1 , dimX*dimT)
    XP <- rep(-1 , dimX*dimT)
    PF <- rep(-1 , dimX*dimX*dimT)
    PP <- rep(-1 , dimX*dimX*dimT)
    YP <- rep(-1 , dimY*dimT)
    R <- rep( -1 , dimY*dimY*dimT)
    KGAIN <- rep(-1 , dimX*dimY*dimT)

    if(echo) {cat("Starting .Fortran() \n")}
    
    # Run compiled code
    FOBJ <-  .Fortran("LTI_KALMAN_FULLA_WITHINPUT",
                      LL   =as.double(LL),
                      INFO =as.integer(INFO),
                      A    =as.double(matA),
                      B    =as.double(matB),
                      C    =as.double(matC),
                      D    =as.double(matD),
                      Time =as.double(Time),
                      Y    =as.double(Y),
                      U    =as.double(U),
                      SIG  =as.double(SIG),
                      S    =as.double(S),
                      X0   =as.double(InitialState),
                      dimT =as.integer(dimT) ,
                      dimY =as.integer(dimY),
                      dimU =as.integer(dimU),
                      dimX =as.integer(dimX),
                      DoseN=as.integer(DoseN),
                      DoseTime=as.double(DoseTime),
                      DoseAmt=as.double(DoseAmt),
                      DoseState=as.integer(DoseState),
                      Xf   =as.double(XF),
                      Xp   =as.double(XP),
                      Pf   =as.double(PF),
                      Pp   =as.double(PP),
                      Yp   =as.double(YP),
                      R    =as.double(R),
                      Kgain=as.double(KGAIN),
                      PACKAGE="PSM")
    
    
    if(echo) cat("\t -iLL: " ,  round(FOBJ$LL,3) , "\n")
    if(outputInternals) {
      R <- array(FOBJ$R,c(dimY,dimY,dimT))
      Yp <- array(FOBJ$Yp,c(dimY,dimT))
      KfGain <- array(FOBJ$Kgain,c(dimX,dimY,dimT))
      R[R>1E200] <- NA
      Yp[Yp>1E200] <- NA
      KfGain[KfGain>1E200] <- NA   
      return(list( negLogLike=FOBJ$LL,
                  Time=array(FOBJ$Time,c(dimT)), 
                  Xp=array(FOBJ$Xp,c(dimX,dimT)),
                  Xf=array(FOBJ$Xf,c(dimX,dimT)),
                  Yp=Yp,
                  KfGain=KfGain,
                  Pf=array(FOBJ$Pf,c(dimX,dimX,dimT)),
                  Pp=array(FOBJ$Pp,c(dimX,dimX,dimT)),
                  R = R
                  )
             )
    } else {
      return(as.vector(FOBJ$LL))
    }
    
    
  } # Not Fast - Ie. singular or user chose R implementation

  
    # P0 CTSM MathGuide page 19 (1.118) and page 8 (1.49)
    PS  <- 1.0;
    tau <- Time[2] - Time[1]

    # Dimensions Pintegral:[2*dimX 2*dimX]
    tmp <-  rbind( cbind(-matA , SIG%*%t.default(SIG)) ,
                  cbind( matrix(rep(0,2*dimX),nrow=dimX,ncol=dimX) , t.default(matA) ))*tau
#print(tmp)
  
    # Use Matrix package to compute Matrix exponential
    Pint  <- matexp(tmp)


  
    PHI0   <- t.default(Pint[(dimX+1):(2*dimX),(dimX+1):(2*dimX),drop=FALSE])
    P0    <- PS*PHI0 %*% Pint[1:dimX,(dimX+1):(2*dimX),drop=FALSE]

    #----------------------------------------------------------
    # Init matrices used in the Kalman filtering
    #----------------------------------------------------------
    Xp      <- array(NA,c(dimX,dimT))
    Xf      <- array(NA,c(dimX,dimT))
    Yp      <- array(NA,c(dimY,dimT))
    KfGain  <- array(NA,c(dimX,dimY,dimT))
    Pf      <- array(NA,c(dimX,dimX,dimT))
    Pp      <- array(NA,c(dimX,dimX,dimT))
    R       <- array(NA,c(dimY,dimY,dimT))

    # Insert initial estimates into matrices
    Pp[,,1] <- P0
    Xp[,1]  <- InitialState



    #----------------------------------------------------------
    # Kalman Filtering for Linear Time-Invariant Model
    #----------------------------------------------------------

    # Init the Negative LogLikelihood
    negLogLike <- 0

    # Variables in the for loop
  DIAGDIMY  <- diag(1,dimY)
  DIAGDIMX  <- diag(1,dimX)
  DIAGRANKA <- diag(1,rankA)
  DIAGDIMXRANKA <- diag(1,dimX-rankA)
  SIGSIGT <- SIG%*%t.default(SIG)
  ZERODIMXDIMX <- matrix(0,nrow=dimX,ncol=dimX)
  matA.T <- t.default(matA)
  matC.T <- t.default(matC)
  if(singA){                                        # matA singular
    Ua <- svd(matA)$u
    Ua.T <- t.default(Ua)    
  } else {
    matA.Inv <- solve.default(matA)
  }


  
  matASIGSIGTzerosmatAT <-  rbind( cbind(-matA , SIGSIGT ) ,
                  cbind( ZERODIMXDIMX , matA.T  ))

    # Loop over timepoints
    for(k in 1:dimT) {
    
        if(echo) cat("Looping \t k= " ,k, ", Time = " , Time[k] , "\n")


        if(echo) cat("Output Predictions...\n")          
        # Does the observation has missing values?
        ObsIndex  <- which(!is.na(Y[,k]))
        E         <- DIAGDIMY[ObsIndex,,drop=FALSE]

        # Set Uk
        if(ModelHasInput) {
          Uk <- U[,k,drop=FALSE]
        } else {Uk <- NA}                

        # Create output prediction
        Yp[ObsIndex,k] <- {
                  if(ModelHasInput) {
                  E %*% (matC%*%Xp[,k,drop=FALSE] + matD%*%Uk)
                  } else { 
                  E%*%matC%*%Xp[,k,drop=FALSE]} }

        if(echo) cat(paste(Yp[ObsIndex,k],"\n",sep="")) 

        # Output prediction covariance    
        #S     <- Model$S(phi=phi)
        R[ObsIndex,ObsIndex,k]  <- E%*%matC%*%Pp[,,k]%*% matC.T %*%t.default(E) + E%*%S%*%t.default(E)

  
        if(length(ObsIndex)>0) { #there must be at least one obs for updating.
        
          if(echo) cat("Updating States ...\n")          
            
          # Kalman gain
          KfGain[,ObsIndex,k] <- Pp[,,k]%*% matC.T %*% t.default(E) %*% solve.default(R[ObsIndex,ObsIndex,k])

          # Updating
          e       <- Y[ObsIndex,k,drop=FALSE]-Yp[ObsIndex,k,drop=FALSE]
          KFg     <- CutThirdDim(KfGain[,ObsIndex,k,drop=FALSE])
          Xf[,k]  <- Xp[,k,drop=FALSE] + KFg%*%e
          Pf[,,k] <- Pp[,,k] - KFg %*% R[ObsIndex,ObsIndex,k] %*% t.default(KFg)
  
          # Add contribution to negLogLike: CTSM page 3. (1.14)
          tmpR        <-  CutThirdDim(R[ObsIndex,ObsIndex,k,drop=FALSE])
          logdetR     <-  determinant.matrix(tmpR)$modulus
     
          negLogLike  <- negLogLike + .5*(logdetR + t.default(e)%*%solve.default(tmpR)%*%e)
        } else { #no observations, update not available.
          Xf[,k] <- Xp[,k]
          Pf[,,k] <- Pp[,,k]
        }

        
        if(DataHasDose) {
          # Check if dosing is occuring at this timepoint.
          if( any(Time[k]==Data$Dose$Time)) {
            if(echo) cat("Dosing...\n")
            idxD = which(Time[k]==Data$Dose$Time)
            # Multiple dosing a timepoint[k]
            for(cmt in 1:length(idxD)) {
              Xf[Data$Dose$State[idxD[cmt]],k] <- Xf[Data$Dose$State[idxD[cmt]],k] + Data$Dose$Amount[idxD[cmt]]
            }
          }
        }
        
        # Abort if negLogLike is Inf or -Inf
        if(is.infinite(negLogLike)) {
          if(echo) cat("\t -LL: " ,  round(negLogLike,3), "\n")
          return(negLogLike)
          }
        
        # Abort state pred if finished
        if(k==dimT) break
        
        if(echo) cat("State Prediction...\n")            

        # State prediction
        tau   <- Time[k+1]-Time[k]

        # Use large time invariant matrix
        if(echo) cat("Starting matrix exponential...\n")
        if(echo) print(tau)
        if(echo) print(matASIGSIGTzerosmatAT)         
        tmp   <- matexp(tau * matASIGSIGTzerosmatAT)
        if(echo) cat("Finished matexp()\n")         
        

        # CTSM (1.48)
        PHI   <- t.default(tmp[(dimX+1):(2*dimX),(dimX+1):(2*dimX),drop=FALSE])
        # CTSM (1.49)
        IntExpASIG <- PHI %*% tmp[1:dimX,(dimX+1):(dimX*2),drop=FALSE]
        # CTSM (1.45)
        Pp[,,k+1] <- PHI %*% Pf[,,k] %*% t.default(PHI) + IntExpASIG

        # Different formulaes depending on A and INPUTS
        if( !singA ) {
            # Special case #3: Non-singular A, zero order hold on inputs.
            Xp[,k+1] <- { if( ModelHasInput) {
                              PHI%*%Xf[,k,drop=FALSE]+ matA.Inv %*%(PHI-DIAGDIMX)%*%matB%*%Uk
                              } else {
                              PHI%*%Xf[,k,drop=FALSE] } }
        } else {

              # Special case #1: Singular A, zero order hold on inputs.
              if(!ModelHasInput) {
                  Xp[,k+1] <- PHI %*% Xf[,k]
                }
              else  {
                  # matA is singular and has INPUT
                  # CTSM Mathguide Special Case no.1 page 10
               
                  PHITilde      <- Ua.T %*% PHI %*% Ua
                  PHITilde1     <- PHITilde[1:rankA,1:rankA,drop=FALSE]
                  # PHITilde2     <- PHITilde[1:rankA,(rankA+1):dimX,drop=F]
                  # PHITilde1Inv  <- solve(PHITilde1)

                  ATilde    <- Ua.T %*% matA %*% Ua
                  ATilde1   <- ATilde[1:rankA,1:rankA,drop=FALSE]
                  ATilde2   <- ATilde[1:rankA,(rankA+1):dimX,drop=FALSE]
                  ATilde1Inv <- solve.default(ATilde1)

                  IntExpAtildeS <- ZERODIMXDIMX
                  # Insert upper left part of matrix [1:rankA 1:rankA]
                  IntExpAtildeS[1:rankA , 1:rankA] <- ATilde1Inv %*% (PHITilde1-DIAGRANKA)
                  # Lower left part
                  # IntExpAtildeS[(rankA+1):dimX , 1:rankA]   <- 0
                  # Upper Right
                  IntExpAtildeS[1:rankA,(rankA+1):dimX] <- ATilde1Inv %*%
                      (IntExpAtildeS[1:rankA,1:rankA]-DIAGRANKA*tau)%*%ATilde2
                  # Lower right
                  # Changed 2007-12-03 SKLI IntExpAtildeS[(rankA+1):dimX,(rankA+1):dimX] <- diag(1,dimX-rankA)*tau
                  IntExpAtildeS[(rankA+1):dimX,(rankA+1):dimX] <- DIAGDIMXRANKA*tau

                  # Insert State prediction CTSM  (1.60)
                  Xp[,k+1] <- PHI %*% Xf[,k]+ Ua %*% IntExpAtildeS %*% Ua.T %*% matB %*% Uk

                  } # end else
            }    # end else
        } #end for

        # Complete negLogLike  with second part CTSM page 3 (1.14)
        #negLogLike  <- negLogLike + .5*dimT*dimY*log(2*pi)
        negLogLike  <- negLogLike + .5*( dimT*dimY - sum(is.na(Y)) )*log(2*pi)

        if(echo) cat("\t -LL: " ,  round(negLogLike,0) , "\n")
        if(outputInternals) {
          return( list( negLogLike=negLogLike,Time=Time,Xp=Xp, Xf=Xf, Yp=Yp, KfGain=KfGain,Pf=Pf, Pp=Pp, R=R))
        } else {
          return(as.vector(negLogLike)) }
}





