`PSM.smooth` <-
function(Model,Data,THETA,subsample=0,trace=0,etaList=NULL) {

  dimS <- length(Data)
  for(i in 1:dimS) {
    check <- ModelCheck(Model,Data[[i]],list(Init=THETA))
    if(!check$ok) {
      errmsg <- check$errmsg
      errmsg <- paste(errmsg, "- the error occured using data for individual", i)
      break
    }
  }
  if(!check$ok) stop(errmsg)
  Linear <- check$Linear
  
  OMEGA <- Model$ModelPar(THETA)$OMEGA
  theta <- Model$ModelPar(THETA)$theta

  FullOutput <- FALSE
  if( is.null(etaList) && !is.null(OMEGA)) {
    apl <- APL.KF(THETA=THETA,Model=Model,Pop.Data=Data,GUIFlag=trace,
                   longOutput=TRUE,Linear=Linear)
    etaList <- apl$etaList
    FullOutput <- TRUE
  }
  
  dimY <- dim(Data[[1]]$Y)[1]
  lkf = vector(mode="list",length=dimS)
  
  for(i in 1:dimS) {
    Di <- Data[[i]]
    Y <- Di$Y
    Time <- Di$Time
    # Check for INPUT and set it
    if(is.null(Di[["U"]])) { #check if U exists.
      ModelHasInput <- FALSE
    } else if ( sum(is.na(Di[["U"]])) >=1 ) { #check if it contains any NA
      ModelHasInput <- FALSE
    } else {
      ModelHasInput <- TRUE
    }
    U <- if( !ModelHasInput) { NA } else { Di[["U"]] }
    if(subsample) {
      n <- length(Time)
      N <- n+(n-1)*subsample
      TT <- vector(length=N)
      YY <- matrix(0,nrow=dimY,ncol=N)

      if(ModelHasInput) {
        UU <- matrix(0,nrow=dim(U)[1],ncol=N)
      } else {
        UU <- NA
      }
      for(j in 1:n) {
        idx0 <- j+(j-1)*subsample
        YY[,idx0] <- Y[,j]
        TT[idx0] <- Time[j]
        if(ModelHasInput)
          UU[,idx0] <- U[,j]
        if(j==n)
          break
        YY[,(idx0+1):(idx0+subsample)] = NA
        TT[(idx0+1):(idx0+subsample)] <- seq(from=Time[j],to=Time[j+1],
                                             length.out=(subsample+2))[-c(1,subsample+2)]
        if(ModelHasInput)
          UU[,(idx0+1):(idx0+subsample)] = matrix(U[,j],nrow=dim(U)[1],ncol = subsample)
      }
      Di <- list(Y = YY, Time = TT, U = UU)
    }
    if(!is.null(OMEGA)) {
      phi <- Model$h(etaList[,i],theta=theta,covar=Data[[i]]$covar)
    } else {
      # OMEGA IS NULL
      phi <- Model$h( eta=NULL , theta=theta , covar=Data[[1]]$covar)
    }
    Di$Dose <- Data[[i]]$Dose
    
    if(Linear) {
      lkf[[i]] <- LinKalmanSmoother( phi=phi, Model=Model , Data=Di )
    } else {
      lkf[[i]] <- ExtKalmanSmoother( phi=phi, Model=Model , Data=Di )
    }
    
    if(trace)
      print(paste("Individual",i))
  

    if(FullOutput) {
      lkf[[i]]$eta <- etaList[,i]
      lkf[[i]]$negLogL <- apl$negLogLike
    } 
  } 

  lkf
}

