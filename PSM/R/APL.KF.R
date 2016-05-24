`APL.KF` <-
function (THETA,Model,Pop.Data,LB=NULL,UB=NULL,GUIFlag=0,longOutput=FALSE,fast=TRUE,Linear=NULL) {
### NOTES -  requires: Model$ModelPar(), 

    # Check argument Linear is set otherwise set it.
  if( is.null(Linear) ) {
    if( "Matrices" %in% names(Model) | "Functions" %in% names(Model)) {
      Linear <- TRUE
      if("Functions" %in% names(Model)) Linear <- FALSE    
    } else {
      stop("Cannot determine if model is Linear or Non-linear based on elements in model. (Matrices or Functions)")
    }
  }

  
  if (!is.null(LB)) {
    THETA <- invlogit(THETA,LB,UB)
  }

  OMEGA <- Model$ModelPar(THETA)$OMEGA
  theta <- Model$ModelPar(THETA)$theta


  # Dimensions
  dimS <- length(Pop.Data)
  ifelse(is.null(OMEGA),dimEta <- 1, dimEta <- dim(OMEGA)[1])

  # Init
  etaList   <- matrix(0,nrow=dimEta,ncol=dimS)
  optimStat <- matrix(0,nrow=3,ncol=dimS)
  LiPart    <- matrix(0,nrow=1,ncol=dimS)
  
  if(GUIFlag>1) 
    starttime <- proc.time()[3]
  

  for (i in 1:dimS) {
    if(GUIFlag>2)
      print(paste('Individual', i))
    if(!is.null(OMEGA)) {
      # OMEGA has a value
      result <- APL.KF.individualloop(theta=theta,OMEGA=OMEGA,Model=Model,
                                      Data=Pop.Data[[i]],GUIFlag=GUIFlag,fast=fast,Linear) 
      LiPart[i] <- result$LiPart_i
      etaList[,i] <- result$eta_i
      optimStat[,i] <- result$optimStat_i
    }
    else {
      # OMEGA is NULL      
      phi <- Model$h(eta=NULL,theta=theta,covar=Pop.Data[[i]]$covar) 

      # Run Kalman Filter according to complexity    
      if(Linear) {
        LiPart[i] <- LinKalmanFilter( phi=phi, Model=Model , Data=Pop.Data[[i]],fast=fast )
      } else {
        LiPart[i] <- ExtKalmanFilter( phi=phi, Model=Model , Data=Pop.Data[[i]] )
      }
      etaList[,i] <- NaN
      optimStat[,i] <- NaN
    }      
  }

  # Boundary Penalty function
  pen <- 0
  if (!is.null(LB)) {
     # Penalty-function
      lambda = 1e-4 #værdi fra ctsm-userguide p. xx
      dx = 1e-30    #
      for(i in 1:length(THETA)) pen = pen + lambda*(abs(LB[i])/(THETA[i]-LB[i]+dx) + abs(UB[i])/(UB[i]-THETA[i]+dx))
  }

  
  if(GUIFlag>1) {
    totaltime = proc.time()[3]-starttime
    minutes = floor(totaltime/60)
    tid <- paste("  (",minutes,":",round(totaltime-60*minutes,2),")",sep="")
    print(c(" -logL  =", signif(sum(LiPart)+pen,10),tid),q=FALSE)

  }
  
  if(longOutput) {
    list(negLogLike=sum(LiPart)+pen,etaList=etaList,optimStat=optimStat)
  } else {  
    #The return variable - neg. Log. Likelihood
    sum(LiPart)+pen
  }

}

