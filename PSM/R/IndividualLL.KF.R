`IndividualLL.KF` <-
function (eta,theta,OMEGA,Model,Data,fast=TRUE,Linear=NULL){
### NOTES -  requires: o$negLogLike, o$Yp, h(eta,theta)

    # Check argument Linear is set otherwise set it.
  if( is.null(Linear) ) {
    if( "Matrices" %in% names(Model) | "Functions" %in% names(Model)) {
      Linear <- TRUE
      if("Functions" %in% names(Model)) Linear <- FALSE    
    } else {
      stop("Cannot determine if model is Linear or Non-linear based on elements in model. (Matrices or Functions)")
    }
  }
  
  phi <- Model$h(eta,theta,covar=Data$covar)

  # run the KF one time, to evaluate negative log-likelihood
  if(Linear) {
    negLogLike <- LinKalmanFilter( phi=phi, Model=Model , Data=Data , fast=fast)
  } else {
    negLogLike <- ExtKalmanFilter( phi=phi, Model=Model , Data=Data )
  } 
  eta <- matrix(eta,ncol=1)

  #Return a posteriori negative log likelihood
  negLogLike + 0.5*t(eta)%*%solve(OMEGA)%*%eta + 0.5*log( det(2*pi*OMEGA))

}

