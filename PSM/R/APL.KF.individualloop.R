`APL.KF.individualloop` <-
function (theta,OMEGA,Model,Data,GUIFlag=0,fast=TRUE,Linear=NULL) {
### NOTES -  requires: o$R, o$Yp, h(eta,theta)


    # Check argument Linear is set otherwise set it.
  if( is.null(Linear) ) {
    if( "Matrices" %in% names(Model) | "Functions" %in% names(Model)) {
      Linear <- TRUE
      if("Functions" %in% names(Model)) Linear <- FALSE    
    } else {
      stop("Cannot determine if model is Linear or Non-linear based on elements in model. (Matrices or Functions)")
    }
  }

  dimEta <- dim(OMEGA)[1]
  tmp <- dim(Data$Y)
  dimY <- tmp[1]
  dimN <- tmp[2]
  eGrad <- array(0,dim = c(dimY,dimEta,dimN))
  
  controllist <- list(trace=0,maxit=500)
  
  if(GUIFlag>2) {
    controllist$trace <- GUIFlag-2  #higher number -> more detail
    controllist$REPORT <- 1 #report for every X iteration
  }
  if(0) { #unconstrained
    controllist$reltol <- 1e-7
    out <- optim(par = rep(0,dimEta), fn = IndividualLL.KF , gr = IndividualLL.KF.gr,
                 method = "BFGS", control = controllist, hessian = FALSE, 
                 theta=theta, OMEGA=OMEGA, Model=Model, Data=Data, fast=fast)
  } else { #constrained
    controllist$factr <- 1e8
    out <- optim(par = rep(0,dimEta), fn = IndividualLL.KF, gr = IndividualLL.KF.gr,
                 method = "L-BFGS-B", lower = -4*sqrt(diag(OMEGA)),upper = 4*sqrt(diag(OMEGA)),
                 control = controllist, hessian = FALSE, 
                 theta=theta, OMEGA=OMEGA, Model=Model, Data=Data, fast=fast, Linear=Linear)
  }
  
  
  # Print optimization stats
  if(GUIFlag>2) {
    print(out$counts)
  }
  optimStat_i <- c(out$convergence, out$value, out$counts[1])
  
  phi <- Model$h(eta=out$par,theta=theta,covar=Data$covar)
  if(Linear) {
    o <- LinKalmanFilter(phi=phi, Model=Model, Data=Data, outputInternals=TRUE, fast=fast)
  } else {
    o <- ExtKalmanFilter(phi=phi, Model=Model, Data=Data, outputInternals=TRUE)
  }
  
  # Calculate stepsize in central diffenrence gradient
  stepSize <- 1E-5;
  
  # Create the e-grad
  for (p in 1:dimEta) {
    d <- rep(0,dimEta)
    d[p] <- stepSize
    
    # Forward difference
    phi.f <- Model$h(out$par+d,theta,covar=Data$covar)
    # Backward difference
    phi.b <- Model$h(out$par-d,theta,covar=Data$covar)

    if(Linear) {
      eF <- Data$Y - LinKalmanFilter(phi=phi.f, Model=Model, Data=Data, outputInternals=TRUE, fast=fast)$Yp
      eB <- Data$Y - LinKalmanFilter(phi=phi.b, Model=Model, Data=Data, outputInternals=TRUE, fast=fast)$Yp
    } else {
      eF <- Data$Y - ExtKalmanFilter(phi=phi.f, Model=Model, Data=Data, outputInternals=TRUE)$Yp
      eB <- Data$Y - ExtKalmanFilter(phi=phi.b, Model=Model, Data=Data, outputInternals=TRUE)$Yp
    }

    #Insert the calculated gradient.
    eGrad[,p,] <- (eF-eB)/(2*stepSize)
  }

  # Hessian approximation (19)
  h_Li <- matrix(0,dimEta,dimEta)
  for (q in 1:dimN) {
    ObsIndex  <- which(!is.na(Data$Y[,q]))
    if(length(ObsIndex)>0) {
      eG <- CutThirdDim(eGrad[ObsIndex,,q,drop=FALSE])
      h_Li <- h_Li + t.default(eG) %*%  solve(o$R[ObsIndex,ObsIndex,q]) %*% eG
    }
  }
  h_Li <- - h_Li - solve(OMEGA);

  # RETURN neg. log. likelihood contribution
  list(   LiPart_i = .5*log(det(-1*h_Li/(2*pi))) + out$value,
             eta_i = out$par,
       optimStat_i = optimStat_i
       )
}

