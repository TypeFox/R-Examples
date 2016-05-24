"kifeltent" <-
function (QKEP=.QKEP,QND=.QND,MIT=c("Z"==substring(names(.QND),1,1))) {

  # INDA IS A VECTOR OF MAXIMUM LEVELS IN EACH OF THE Z VARIABLES

  cat("\nEntering kifeltent()")
  
  if(length(QND)!=length(MIT)) {
    return(NULL)
  }
  INDA <- QND[MIT]
  HI <- length(INDA)
  IND <- NULL

  cat("\nINDA: ",INDA,"\nHI: ",HI, "\nIND: ",IND)

  for (hi in 1:HI){
    IND <- desca((1:INDA[hi]),IND)
  }
  
  N <- dim(IND)[1]
  if(is.null(N)) {
    N <- length(IND)
  
    if(is.null(N)) {
      return()
    }
    else {
      IND<-matrix(IND,ncol=1)
    }
  }

  KI <- matrix(nrow=N,ncol=2,dimnames=list(NULL,c("HCONDI","PCONDI")))

  for (n in 1:N) {
    Q <- szelet(QKEP,MIT,IND[n,])
    SQ <- sum(Q)
    if(SQ > 0) {
      KI[n,1] <- entro(Q/SQ)
      KI[n,2] <- SQ
    }
    else {
      KI[n,1] <- 0
      KI[n,2] <- 0
    }

  }
  
  KI
   
  cat("\nLeaving kifeltent()")
  
}

