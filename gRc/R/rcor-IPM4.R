## IPM for RCOR models
##
## Fit model where K=ACA. Hence logL is of the form
## logL = |K| - tr(KS) = |K| - tr(ACAS)
## 1. Start with current estimate of C, e.g. C=I
## 2. For known C maximize logL over paramters in A.
##    Is done by iteratively solving system of 2nd order equations using refitA
## 3. For known A, write tr(ACAS)=tr(CASA)=tr(CQ), where Q=ASA
##    Is done by regarding Q=ASA as "data matrix" and then estimating elements
##    C using IPM only on off-diagonals of C
##


rcorIPM <- function(object, K0, control=object$control, trace=object$trace){
  
  S     <- object$dataRep$S
  n     <- object$dataRep$n
  VCC   <- object$intRep$vccI
  glist <- object$intRep$eccI

  logL0       <- prevlogL <- ellK(K0,S,n-1)
  
  maxit       <- control$maxouter
  #cat("maxouter:", maxit,"\n")
  logLeps     <- control$logLeps * abs(logL0)
  #cat("logLeps:", logLeps, "\n")
  logL.vec    <- rep(NA, maxit)
  converged   <- FALSE

  ctrl      <- control
  deltaeps  = ctrl$deltaeps
  maxouter  = ctrl$maxouter
  maxinner  = ctrl$maxinner
  nobs      = n

  logLINT      = 0
  convergedINT = 1

  C.curr    <- cov2cor(K0)             ## print(C.curr)
  A         <- refitA(S, C.curr, VCC);   ## First A estimate

                                        #cat("maxouter:", maxouter,"\n")
  itcount <- 1
  while(!converged){
    Q.curr        <- A * t((A*S)) ## Short for diag(A) %*% S %*% diag(A)
    
    ##cat("rcor - call C\n")
    ansC <- .Call("rconipm", S=Q.curr, nobs=nobs-1, C.curr, Glist=glist, 
                  maxouter=maxouter, maxinner=maxinner, 
                  logL=logLINT, logLeps=logLeps, deltaeps=deltaeps,
                  converged=convergedINT, debug=trace,
                  PACKAGE="gRc")

    C.curr<- ansC[[1]]

    K.curr  <- A * t((A*C.curr)) ## Short for diag(A) %*% C.curr %*% diag(A)

    logL    <- ellK(K.curr,S,n-1);
    dlogL   <- logL - prevlogL
    
    if (trace>=2)
      cat("..rcorIPM iteration (offdiag) ", itcount, "logL:", logL, "dlogL:", dlogL, "\n")
    
    A       <- refitA(S, C.curr, VCC, Astart=A)
    ##     KKKK  <- A * t((A*C.curr)) ## Short for diag(A) %*% C.curr %*% diag(A)
    ##     logLLLL    <- ellK(KKKK,S,n-1);
    ##     if (trace>=2)
    ##       cat("..rcorIPM iteration (diag)    ", itcount, "logL:", logLLLL, "dlogL:", dlogL, "\n")
    
    logL.vec[itcount]  <- logL

    if (dlogL < logLeps || itcount>=maxit)
      converged <- TRUE
    else {
      prevlogL <- logL;
      itcount  <- itcount + 1
    }
  }
  
  coef <- K2theta(object,K.curr, scale='original')
  vn   <- unlist(lapply(getcc(object),names))
  names(coef) <- vn

  if (!is.null(control$vcov)){
    J  <- getScore(object,K=K.curr)$J
    dimnames(J) <- list(vn, vn)
  } else {
    J <- NULL
  }
  
  dimnames(K.curr) <- dimnames(S)
  ans <- list(K=K.curr, logL=logL, coef=coef, J=J, logL.vec=logL.vec[1:itcount])
  return(ans)
}


refitA <- function(S, K, vccI, Astart=NULL, itmax=100){
  
  d <- nrow(S)
  if (is.null(Astart))
    Astart <- rep(1,d)
  Aprev <- Astart
  iii   <- 0
  
  all <- unique(unlist(vccI))
  all <- all[order(all)]

  gset <- vccI
  cset <- vccI

  for (i in 1:length(vccI)){
    gset[[i]] <- as.numeric(vccI[[i]])
    cset[[i]] <- setdiffPrim(all, gset[[i]])
  }

  repeat{
    iii <- iii + 1
    for (i in 1:length(vccI)){
      cns   <- gset[[i]]
      compl <- cset[[i]]
      Ai <- sum(S[cns,cns] * K[cns,cns])
      Bi <- sum(S[compl,cns] * K[compl,cns] *Astart[compl])
      Ci <- length(cns)
      Di <- (Bi^2 + 4*Ai*Ci)
      xi <- (-Bi + sqrt(Di))/(2*Ai)
      Astart[cns] <- xi
      ##cat("Ai:", Ai, "Bi:", Bi, "Ci:", Ci, "\n")
    }
    if ((max( (Aprev-Astart)^2) < 1e-4) | iii>=itmax){
      ##cat("Iterations:",iii,"\n\n")
      break
    }    
    Aprev <- Astart
  }
  attr(Astart,"iterations") <- iii
  return(Astart)
}



# ## OLD ONE
# .rcorIPM <- function(object, K0, control=object$control, trace=object$trace){
#   S   <- dataRep(object, "S")
#   n   <- dataRep(object, "n")
#   VCC <- intRep(object,"vccI")
#   V   <- NULL

#   logL0       <- prevlogL <- ellK(K0,S,n-1)
  
#   maxit       <- control$maxouter
#   logLeps     <- control$logLeps * abs(logL0)

#   logL.vec    <- rep(NA, maxit)
#   converged   <- FALSE

#   ## Temporary model object with new controls
#   ctrl        <- object$control
#   ctrl$vccfit <- FALSE
#   ctrl$vcov   <- NULL
#   o2          <- object
#   o2$control  <- ctrl

#   ##C.curr    <- diag(1, nrow(S));     ## Initial C
#   C.curr    <- cov2cor(K0)             ##print(C.curr)
#   A         <- refitA(S,C.curr,VCC); ## First A estimate

#   itcount <- 1
#   while(!converged){
#     Q.curr        <- A * t((A*S)) ## Short for diag(A) %*% S %*% diag(A)
#     o2$dataRep$S  <- Q.curr

#     C.curr  <- rconIPM(o2,K0=C.curr, control=o2$control, trace=trace-1)$K
#     K.curr  <- A * t((A*C.curr)) ## Short for diag(A) %*% C.curr %*% diag(A)
#     A       <- refitA(S, C.curr, VCC, Astart=A)
#     logL    <- ellK(K.curr,S,n-1);
#     dlogL   <- logL - prevlogL
    
#     if (trace>=3)
#       cat("...rcorIPM iteration", itcount, "logL:", logL, "dlogL:", dlogL, "\n")

#     logL.vec[itcount]  <- logL
#     if ((logL-prevlogL) < logLeps || itcount>=maxit)
#       converged <- TRUE
#     else {
#       prevlogL <- logL;
#       itcount  <- itcount + 1
#     }
#   }

#   coef <- K2theta(object,K.curr, scale='original')
#   vn   <- unlist(lapply(getcc(object),names))
#   names(coef) <- vn

#   if (!is.null(control$vcov)){
#     J  <- getScore(object,K=K.curr)$J
#     dimnames(J) <- list(vn, vn)
#   } else {
#     J <- NULL
#   }
   
#   dimnames(K.curr) <- dimnames(S)
#   ans <- list(K=K.curr, logL=logL, coef=coef, J=J, logL.vec=logL.vec[1:itcount])
#   return(ans)
# }


# refitA <- function(S, K, NSi, Astart=NULL, itmax=100){

#   cstart <- Astart
#   d <- nrow(S)
#   if (is.null(cstart))
#     cstart <- rep(1,d)
#   cprev <- cstart
#   iii   <- 0

#   all <- unique(unlist(NSi))
#   all <- all[order(all)]
#   repeat{
#     iii <- iii + 1
#     for (i in 1:length(NSi)){
#       cns   <- as.numeric(NSi[[i]])

#       compl <- setdiff(all,cns)
#       Ai <- sum(S[cns,cns] * K[cns,cns])
#       Bi <- sum(S[compl,cns] * K[compl,cns] *cstart[compl])
#       Ci <- length(cns)
#       Di <- (Bi^2 + 4*Ai*Ci)
#       xi <- (-Bi + sqrt(Di))/(2*Ai)
#       cstart[cns] <- xi
#       ##cat("Ai:", Ai, "Bi:", Bi, "Ci:", Ci, "\n")
#     }
#     if ((sum( (cprev-cstart)^2) < 1e-9) | iii>=itmax){
#       ##cat("Iterations:",iii,"\n\n")
#       break
#     }
#     cprev <- cstart
#   }
#   attr(cstart,"iterations") <- iii
#   return(cstart)
# }


# rcorFitIterative <- function(m,control=rcox.control()){

#   if (control$trace>=3)
#     cat("...rcorFitIterative\n")
  
#   tstart <- proc.time()
#   conv   <- control$logLepsilon
#   itmax  <- control$maxit
#   rconModel     <- m$call
#   method <- rconModel$method
  
#   VCC       <- c(m$stdrep$VCCU, m$stdrep$VCCR)
#   nam       <- m$varnames
#   n         <- m$n
#   S         <- m$S
#   C.curr    <- diag(1, nrow(S));  #print("C.curr (Start)"); print(C.curr)
#   ll        <- d.ll <- NULL
#   prev.logL <- iii<- 0
#   a         <- rescaleC(S,C.curr,VCC); #print("a"); print(a)
  
#   controlNew                <- control
#   controlNew$VCCU           <- controlNew$VCCR <- FALSE
#   controlNew$representation <- 'stdrep'
#   controlNew$method         <- 'user'

#   rconModel$type    <- 'rcon'
#   rconModel$method  <- 'user'
#   rconModel$data    <- NULL
#   rconModel$n       <- n
#   rconModel$control <- controlNew
#   rconModel$fit     <- TRUE
#   ##print(rconModel)
  
#   repeat{    
#     iii <- iii + 1
   
#     S.curr           <- diag(a) %*% S %*% diag(a)
#     dimnames(S.curr) <- dimnames(S)
#     ##print("S.curr"); print(S.curr)
    
#     ## Estimate rho
#     rconModel$Kstart <- C.curr
#     rconModel$S      <- S.curr    
#     rconModelnew     <- eval(rconModel)
#     C.curr           <- rconModelnew$fit$K;  ##print("C.curr"); print(C.curr)
#     ## Calculate K
#     K.curr    <- diag(a) %*% C.curr %*% diag(a)

#     ## Estimate a
#     a         <- rescaleC(S, C.curr, VCC)
    
#     ## Monitoring convergence
#     curr.logL <- ellK(K.curr, S, n)
#     ll        <- c(ll, curr.logL)
#     diff.logL <- curr.logL-prev.logL
#     prev.logL <- curr.logL
#     d.ll      <- c(d.ll, diff.logL)

#     if (control$trace>=3)
#       cat ("...rcorFit (iteration)! logL:", curr.logL, "\n");

#     if ((abs(diff.logL)/abs(prev.logL) < conv) || iii>=itmax){
#       break
#     }
#   }
#   dimnames(K.curr) <- dimnames(C.curr)
#   usedTime <- (proc.time()-tstart)[3]
#   value <- list(K=K.curr, logL=curr.logL, logLvec=ll, time=usedTime)

#   if (control$trace>=3)
#     cat ("...rcorFit done! logL:", curr.logL, "\n");
#   return(value)
# }



# rescaleC <- function(S,K,NSi,cstart=NULL,itmax=100){
#   #print(S); print(K); print(NSi)
  
#   d <- nrow(S)
#   if (is.null(cstart))
#     cstart <- rep(1,d)
#   cprev <- cstart
#   iii   <- 0

#   all <- unique(unlist(NSi))
#   all <- all[order(all)]
#   repeat{
#     iii <- iii + 1
#     for (i in 1:length(NSi)){
#       cns   <- NSi[[i]]
#       compl <- setdiff(all,cns)
#       Ai <- sum(S[cns,cns] * K[cns,cns])
#       Bi <- sum(S[compl,cns] * K[compl,cns] *cstart[compl])
#       Ci <- length(cns)
#       Di <- (Bi^2 + 4*Ai*Ci)
#       xi <- (-Bi + sqrt(Di))/(2*Ai)
#       cstart[cns] <- xi
#       ##cat("Ai:", Ai, "Bi:", Bi, "Ci:", Ci, "\n")
#     }
#     if ((sum( (cprev-cstart)^2) < 1e-9) | iii>=itmax){
#       #cat("Iterations:",iii,"\n\n")
#       break
#     }
#     cprev <- cstart
#   }
#   attr(cstart,"iterations") <- iii
#   return(cstart)
# }


