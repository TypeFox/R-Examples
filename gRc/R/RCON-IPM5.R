
rconIPM <- function(object, K0,
                     control=object$control, trace=object$trace){

  if (trace>=2)
    cat("..fitting with rconIPM (C version):\n"); 

  t0 <- proc.time()
  S       <- object$dataRep$S
  nobs    <- object$dataRep$n
  Kstart  <- object$Kstart


  ctrl <- control
  logL      = 0
  converged = 1

  deltaeps  = ctrl$deltaeps
  maxouter  = ctrl$maxouter
  maxinner  = ctrl$maxinner
  
  eccfit <- control$eccfit
  vccfit <- control$vccfit

  if (vccfit)
    vcc <- object$intRep$vccI
  else
    vcc <- NULL

  if (eccfit)
    ecc <- object$intRep$eccI
  else
    ecc <- NULL

  glist   <- c(vcc,ecc)

  logL0       <- prevlogL <- ellK(K0,S,nobs-1)
  logLeps     <- ctrl$logLeps * abs(logL0)

  t00 <- proc.time()
  ##cat("maxouter:", maxouter,"\n")
  Kwork<-.Call("rconipm", S=S, nobs=nobs-1, K=K0, Glist=glist, 
               maxouter=maxouter, maxinner=maxinner, 
               logL=logL, logLeps=logLeps, deltaeps=deltaeps,
               converged=converged, debug=trace,
               PACKAGE="gRc")[[1]]
  ##Kwork <- ansC[[1]]
  ##cat("C-call:", proc.time()-t00, "\n")
  #cat("maxouter:", maxouter, "\n")
  
  coef <- K2theta(object,Kwork, scale='original')
  vn   <- unlist(lapply(getcc(object),names))
  names(coef) <- vn

  if (object$method=="ipms"){ ## IPM without finalizing with finding score
    J <- NULL
  } else {
    if (!is.null(control$vcov)){
      J  <- getScore(object,K=Kwork)$J
      dimnames(J) <- list(vn, vn)
    } else {
      J <- NULL
    }
  }

  if (trace>=2)
    cat("..ipm, logL:", logL, "Time:", proc.time()-t0, "\n")
  
  ##cat("exit rconIPM: Kwork: \n"); print(Kwork)
  ans <- list(K=Kwork, logL=logL, coef=coef, J=J)
  return(ans)
}


