
fit.rcox <- function(object,
                     Kstart  = object$Kstart,
                     method  = object$method,
                     control = object$control,
                     details = object$details,
                     trace   = object$trace,
                     returnModel=TRUE,...){

##  cat("fit.rcox\n")
  if (is.null(Kstart)){
    #cat("Finding Kstart\n")
    Kstart    <- matching(object, trace=trace)$K
    #print(Kstart)
  }

  ##cat("Kstart:\n"); print(Kstart)
  ##Kstart    <- findKinModel(object, KS=object$Kstart,type=object$type, regularize=TRUE)

  ##print(method); print(trace)
  tstart <- proc.time()
  ans <- switch(method,
                "matching"=
                {
                  scoring(object, K0=Kstart, control=control, maxit=1, trace=trace)
                },
                "scoring"=,
                "ipm" =,
                "ipms"=
                {
                  switch(method,
                         "scoring"={
                                        #print(Kstart)
                           scoring(object, K0=Kstart, control=control, trace=trace)
                         },
                         "ipm"=,
                         "ipms"={
                                        #print(Kstart)
                           zz<- ipm(object, K0=Kstart, control=control, trace=trace)         
                           zz
                         })
                },
                "hybrid1"={
                  object2 <- object
                  ctrl          <- object$control
                  ctrl$maxouter <- ctrl$hybrid1switch
                  ctrl$vcov     <- NULL
                  KK  <-ipm(object2, K0=Kstart, control=ctrl, trace=trace)$K
                  scoring(object, K0=KK, control=control, trace=trace)
                }
                )
  ans$method <- method
  ans$Kstart <- Kstart 
  ans$time   <- (proc.time()-tstart)[3]
  
  if (returnModel){
    object$Kstart   <- ans$Kstart
    ans$Kstart <- NULL    
    object$fitInfo  <- ans
    object$method   <- method
    return(object)
  } else {
    return(ans)
  }
}

matching      <- function(object, control=object$control, trace=object$trace){
  if (inherits(object,"rcon"))
    rconScoreMatch(object, control=control, trace=trace)
  else
    rcorScoreMatch(object, control=control, trace=trace)
  ##UseMethod("matching")
}

ipm <- function(object, K0, control=object$control, trace=object$trace){
  if (inherits(object,"rcon"))
    rconIPM(object, K0, control, trace)
  else
    rcorIPM(object, K0, control, trace)
  ##UseMethod("ipm")
}


#scoring <- function(object, K0, control=object$control, maxit=control$maxouter,trace=object$trace) {
#  UseMethod("scoring")
#}

matching.rcon <- function(object, control=object$control, trace=object$trace){
  rconScoreMatch(object, control=control, trace=trace)
}

matching.rcor <- function(object, control=object$control, trace=object$trace){
  rcorScoreMatch(object, control=control, trace=trace)
}

ipm.rcon <- function(object, K0, control=object$control, trace=object$trace){
  rconIPM(object, K0, control, trace)
}

ipm.rcor <- function(object, K0, control=object$control, trace=object$trace){
  rcorIPM(object, K0, control, trace)
}









