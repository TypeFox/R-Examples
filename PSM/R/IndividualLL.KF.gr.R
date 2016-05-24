`IndividualLL.KF.gr` <-
function (eta,theta,OMEGA,Model,Data, GradSTEP=1e-4,GUIFlag=0, fast=TRUE, Linear=NULL) {
  # Forward gradient function for IndividualLL.KF

  # Check argument Linear is set otherwise set it.
  if( is.null(Linear) ) {
    if( "Matrices" %in% names(Model) | "Functions" %in% names(Model)) {
      Linear <- TRUE
      if("Functions" %in% names(Model)) Linear <- FALSE    
    } else {
      stop("Cannot determine if model is Linear or Non-linear based on elements in model. (Matrices or Functions)")
    }
  }    
  
  L <- length(eta)
  ILL <- rep(0,length=(L+1))
  GRAD <- rep(0,length=L)
  
  TP <- matrix(0,nrow=L+1,ncol=L)
  for (i in 1:L) {
    TP[i, ] <- eta
    TP[i,i] <- eta[i] + GradSTEP * ( abs(eta[i])  + GradSTEP );
  }
  TP[L+1,] <- eta
  
  for ( i in 1:(L+1) )
    ILL[i] <- IndividualLL.KF(eta=TP[i,],theta=theta,OMEGA=OMEGA,Model=Model,Data=Data,fast=fast,Linear)
        
  # Calculate Gradient and insert
  for (i in 1:L)
    GRAD[i] = (ILL[i]-ILL[L+1])/( TP[i,i]-eta[i] ) 

  if(GUIFlag>4) {
    print(c("<eta> =", signif(eta,5)),q=FALSE)
    print(c("<GR>    =", signif(GRAD,5)),q=FALSE)
  }
   
  GRAD
}

