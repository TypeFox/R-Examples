`PSM.simulate` <- 
function(Model, Data, THETA, deltaTime, longX=TRUE) {

  # Extract Dimenstions of Data and Check
  dimS <- length(Data)
  if(dimS<1) {
    print("Length of Data is less than 1.")
    break
  }
  
  # Check Model and each Data element  
  for(i in 1:dimS) {
    check <- ModelCheck(Model,Data[[i]],list(Init=THETA),DataHasY=FALSE)
    if(!check$ok) {
      errmsg <- check$errmsg
      errmsg <- paste(errmsg, "- the error occured using data for individual", i)
      break
    }
  }
  if(!check$ok) stop(errmsg)

  Linear = check$Linear
  # if(!Linear) stop('Simulation only implemented for linear models.')


  Result <- Tlist <- Ulist <- covarlist <- vector(mode="list",length=dimS)
  for (i in 1:dimS) {
    Tlist[[i]] <- Data[[i]]$Time
    Ulist[[i]] <- Data[[i]]$U
    covarlist[[i]] <- Data[[i]]$covar
  }

  if(is.null(Data[[1]]$U)) { #No input present
    Ulist <- NULL
  }  
  tmp   <- Model$ModelPar(THETA)
  OMEGA <- tmp$OMEGA
  theta <- tmp$theta

  if(!is.null(OMEGA)) {
    dimEta <- dim(OMEGA)[1]
    # New implementation based on MASS mvrnorm
    eta <- t.default(mvrnorm(n=dimS , mu=rep(0,dimEta) , Sigma=OMEGA))
  } else {
    eta <- NULL
  }

  cat("Simulating individual: ")
  for (i in 1:dimS) {
    cat(paste(i,', ',sep=""))
    if(!is.null(OMEGA)) {
      phi <- Model$h(eta=eta[,i],theta=theta,covar=covarlist[[i]])
    }else {
      # OMEGA IS NULL
      phi <- Model$h( eta=NULL , theta=theta , covar=covarlist[[i]])
    }    

    SampleTime <- Tlist[[i]]
    len <- length(SampleTime)
    tseq <- seq(from = SampleTime[1],to=SampleTime[len],by=deltaTime)
    TimeSpan <- SampleTime[len]-SampleTime[1]
    n <- as.integer(round(TimeSpan/deltaTime + 1))
    
    # find where SampleTime are equal to t
    where <- abs(t.default(matrix(rep(SampleTime,n),nrow=len,ncol=n))-tseq) < 1000*.Machine$double.eps
    
    # Check subsampling and choice of deltaTime
    if(sum(where)!=len)
      stop(simpleError("All sample times must belong to t(0)+n*deltaTime, where n is an integer."))
    
    # Check for INPUT and subsample U
    if( is.null(Ulist) ) { #check if U exists.
      ModelHasInput <- FALSE
    } else {
      ModelHasInput <- TRUE
    }
    
    #check if it contains any NA
    if( ModelHasInput) if (any(is.na(Ulist)) ) stop("Input contains NA.")
    
    if( ModelHasInput) {
      RepIdx <- rep( 1:length(SampleTime) , times=round(c(diff(SampleTime)/deltaTime,1)))
      U <- Ulist[[i]][ , RepIdx ,drop=FALSE]
      Ustart <- U[,1,drop=FALSE]
    } else {
      Ustart <- NA
    }

    # Create logical depending on Dose    
    DataHasDose <- "Dose" %in% names(Data[[i]])

    # Initial States
    InitialState  <- Model$X0(Time=t[1], phi=phi, U = Ustart)
    
    # Dimensions
    dimX <- nrow(InitialState)

    if(Linear) {
      # Linear simulation 
      # Create Matrices
      tmpM  <- Model$Matrices(phi=phi)
      matA  <- tmpM$matA
      matB  <- tmpM$matB
      matC  <- tmpM$matC
      matD  <- tmpM$matD 

      dimY <- nrow(matC)

      # Model components
      SIG   <- Model$SIG(phi=phi)  
      Scov  <- Model$S(phi=phi)
    
      # Is A singular ???
      rankA <- qr(matA)$rank
      singA <- (rankA<dimX)
      
      tmp     <- rep(0, max(dimX,dimY)*n)
      X       <- tmp[1:(dimX*n)]
      dim(X)  <- c(dimX,n)
      Y       <- tmp[1:(dimY*n)]
      dim(Y)  <- c(dimY,n)
  
      # Variables in the for loop
      DIAGDIMX  <- diag(1,dimX)
      DIAGRANKA <- diag(1,rankA)
      BigMat    <- rbind(cbind(-matA , SIG%*%t.default(SIG)) ,
                        cbind( matrix(0,nrow=dimX,ncol=dimX) , t.default(matA) ))
      matDIMXDIMX     <- matrix(NA,dimX,dimX)
    

      # Wiener and Gaussian noise - Not Time dependent
      eW <-   matrix(rnorm(n*dimX),nrow=dimX,ncol=n)
      eObs <- t.default(mvrnorm(n=n*dimX , mu=rep(0,dimY), Sigma=Scov))
    
      # Store Initial State
      X[,1] <- InitialState
    
      # Start looping through timepoints   
      for (k in 1:n) {
 
        if(ModelHasInput) Uk <- U[,k,drop=FALSE]

        # observation
        if(ModelHasInput) {
          Y[,k] <- matC %*% X[,k,drop=FALSE] + matD%*%Uk + eObs[,k]
        } else {
          Y[,k] <- matC %*% X[,k,drop=FALSE] + eObs[,k]
        }
      
        # Add dose after measurement is taken at Time[k]
        if(DataHasDose) {
          idxD = which(tseq[k]==Data[[i]]$Dose$Time)
          if(length(idxD)==1) {
            X[Data[[i]]$Dose$State[idxD],k] <- X[Data[[i]]$Dose$State[idxD],k] + Data[[i]]$Dose$Amount[idxD]
          }
        }
      
 
        # Abort state pred if finished
        if (k == n) break

        # State prediction
        tmp <- deltaTime * BigMat
        tmp <- matexp(tmp)
        PHI <- t.default( tmp[(dimX+1):(2*dimX),(dimX+1):(2*dimX),drop=FALSE])
        IntExpASIG <- PHI %*% tmp[1:dimX,(dimX+1):(dimX*2),drop=FALSE]

    # Different formulaes depending on A and INPUTS
      if( !singA ) {
      # Special case #3: Non-singular A, zero order hold on inputs.
        X[,k+1] <- { if( ModelHasInput) {
          PHI%*%X[,k,drop=FALSE] + solve.default(matA)%*%(PHI-DIAGDIMX)%*%matB%*%Uk
        } else {
          PHI%*%X[,k,drop=FALSE] } }
      } else {
      # Special case #1: Singular A, zero order hold on inputs.
        if(!ModelHasInput) {
          X[,k+1] <- PHI %*% X[,k,drop=FALSE]
        }
        else  {
          if(rankA==0) { #only the zero matrix has rank 0, thus A=0: dx=b*u*deltaTime
            X[,k+1] <- X[,k] +  matB*Uk*deltaTime
          } else {
          
                  # matA is singular and has INPUT and dimX>1
                  # CTSM Mathguide Special Case no.1 page 10
            Ua <- svd(matA)$u
            PHITilde      <- t.default(Ua) %*% PHI %*% Ua
            PHITilde1     <- PHITilde[1:rankA,1:rankA,drop=FALSE]
                  # PHITilde2     <- PHITilde[1:rankA,(rankA+1):dimX,drop=F]
                  # PHITilde1Inv  <- solve(PHITilde1)
            
            ATilde    <- t.default(Ua) %*% matA %*% Ua
            ATilde1   <- ATilde[1:rankA,1:rankA,drop=FALSE]
            ATilde2   <- ATilde[1:rankA,(rankA+1):dimX,drop=FALSE]
            ATilde1Inv <- solve.default(ATilde1)
            
            IntExpAtildeS <- matDIMXDIMX
                  # Insert upper left part of matrix [1:rankA 1:rankA]
            IntExpAtildeS[1:rankA , 1:rankA] <- ATilde1Inv %*% (PHITilde1-DIAGRANKA)
                  # Lower left part
            IntExpAtildeS[(rankA+1):dimX , 1:rankA]   <- 0
                  # Upper Right
            IntExpAtildeS[1:rankA,(rankA+1):dimX] <- ATilde1Inv %*%
              (IntExpAtildeS[1:rankA,1:rankA]-DIAGRANKA*deltaTime)%*%ATilde2
                  # Lower right
            IntExpAtildeS[(rankA+1):dimX,(rankA+1):dimX] <- diag(1,dimX-rankA)*deltaTime

                  # Insert State prediction CTSM  (1.60)
            X[,k+1] <- PHI %*% X[,k]+ Ua %*% IntExpAtildeS %*% t.default(Ua) %*% matB %*% Uk
          } # end if (A is zero matrix) else ...
        } # end else
      } # end else
      
    # Wiener noise
      ew <- SIG %*% (sqrt(deltaTime) * eW[,k+1])
      X[,k+1] <- X[,k+1] + ew
    
    } #end for


      # End Linear Simulation
    } else {
      #########################
      # Non-Linear simulation
      #########################

      # Extract functions
      f  <- Model$Functions$f
      df <- Model$Functions$df
      g  <- Model$Functions$g
      dg <- Model$Functions$dg
      
      # Determine Y-dimension      
      dimY <- length( g(x=InitialState,u=Ustart,time=tseq[1],phi=phi) )

      # Initialize Arrays
      tmp     <- rep(0, max(dimX,dimY)*n)
      X       <- tmp[1:(dimX*n)]
      dim(X)  <- c(dimX,n)
      Y       <- tmp[1:(dimY*n)]
      dim(Y)  <- c(dimY,n)
      
      # Wiener Noise noise
      eW <-   matrix(rnorm(n*dimX),nrow=dimX,ncol=n)
    
      # Insert X0 in X-array
      X[,1] <- InitialState
      
      # Start looping through timepoints   
      for (k in 1:n) {
        
        # Update input to present timepoint
        if(ModelHasInput) Uk <- U[,k,drop=FALSE]
        
        # observation
        Y[,k] <- g(x=X[,k,drop=FALSE],u=Uk,time=tseq[k],phi=phi)        
        # Gaussian Noise
        ObsErr <- t.default(mvrnorm(n=1 , mu=rep(0,dimY), Sigma=Model$S(u=Uk,time=tseq[k],phi=phi)))
        # Obs + Noise
        Y[,k] <- Y[,k] + ObsErr
        
        # Add dose after measurement is taken at tseq[k]
        if(DataHasDose) {
           # Check if dosing is occuring at this timepoint.
          if( any(tseq[k]==Data[[i]]$Dose$Time)) {
            idxD = which(tseq[k]==Data[[i]]$Dose$Time)
            # Multiple dosing a timepoint[k]
            for(cmt in 1:length(idxD)) {
              X[Data[[i]]$Dose$State[idxD[cmt]],k] <-
                X[Data[[i]]$Dose$State[idxD[cmt]],k] + Data[[i]]$Dose$Amount[idxD[cmt]]
            }
          }
        }
        
        # Abort state pred if finished
        if (k == n) break
        
        # State Prediction
        X[,k+1] <- X[,k] + deltaTime*f(x=X[,k],u=Uk,time=tseq[k],phi=phi)
        
        # Add Wiener noise
        syserr <- Model$SIG(u=Uk,time=tseq[k],phi=phi) %*% (sqrt(deltaTime) * eW[,k+1])
        X[,k+1] <- X[,k+1] + syserr

      } # End Looping over timepoints
      
    }  # End Non-Linear simulation
    
    # SubSampling
    idx <- which(( where %*% rep(TRUE,len)) == 1)

    Result[[i]] <- list(X=matrix(X[,idx],nrow=dimX),Y=matrix(Y[,idx],nrow=dimY),
                        Time=SampleTime,U=Ulist[[i]],eta=eta[,i])
    if(DataHasDose) {
      Result[[i]]$Dose <- Data[[i]]$Dose
    }
    if(longX) {
      Result[[i]]$longX <- X
      Result[[i]]$longTime <- tseq
    }
  } #end individual loop
  cat("Done\n")
  # list(Xlist=Xlist,Ylist=Ylist,Tlist=SampleTlist,Ulist=Ulist,eta=eta)
  return(Result)
}
