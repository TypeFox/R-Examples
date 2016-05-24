`APL.KF.gr` <-
function (THETA,Model,Pop.Data,LB=NULL,UB=NULL,GradSTEP=1e-4,GUIFlag=0,fast=TRUE,Linear=NULL) {
  # Forward gradient function for APL.KF

    # Check argument Linear is set otherwise set it.
  if( is.null(Linear) ) {
    if( "Matrices" %in% names(Model) | "Functions" %in% names(Model)) {
      Linear <- TRUE
      if("Functions" %in% names(Model)) Linear <- FALSE    
    } else {
      stop("Cannot determine if model is Linear or Non-linear based on elements in model. (Matrices or Functions)")
    }
  }
  
  L <- length(THETA)
  APL <- rep(0,length=(L+1))
  GRAD <- rep(0,length=L)
  
  TP <- matrix(0,nrow=L+1,ncol=L)
  for (i in 1:L) {
    TP[i, ] <- THETA
    TP[i,i] <- THETA[i] + GradSTEP * ( abs(THETA[i])  + GradSTEP );
  }
  TP[L+1,] <- THETA
  
  for ( i in 1:(L+1) )
    APL[i] <- APL.KF(TP[i,],Model=Model,Pop.Data=Pop.Data,LB=LB,UB=UB,GUIFlag=0,fast=fast,Linear=Linear)
        
  # Calculate Gradient and insert
  for (i in 1:L)
    GRAD[i] = (APL[i]-APL[L+1])/( TP[i,i]-THETA[i] ) 

  if(GUIFlag>1) {
    if(!is.null(LB))
      print(c(" THETA  =", signif(invlogit(THETA,LB,UB),5)),q=FALSE)
    print(c("<THETA> =", signif(THETA,5)),q=FALSE)
    print(c("<GR>    =", signif(GRAD,5)),q=FALSE)
  }
   
  GRAD
}

