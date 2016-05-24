## Multivariate Newton
## 

## Start
#   Get start value for theta(0) = fn(K(0))
#   Get f
#   Calc: logL(0)
#

## At iteration p; given theta(p-1) and K(p-1)
#   Set:  stepsize = 1
#   Get (score, information) = f(K,...)
#   Calc: updating 'discrepancy' DISC 

## Linesearch:
#   Calc: theta(p) = theta(p-1) + DISC*stepsize
#   Calc: K(p) = fn(theta(p))
#   Calc: logL(p)

# If logL(p)<logL(p-1)
#   Set: stepsize = stepsize / 2
#   Do linesearch

# Stop doing lineserach when steps have been halved more than 3 times...

##scoring.rcox <- function(m, K0,
                         
scoring <- function(object, K0, 
                    control = object$control,
                    maxit   = control$maxouter,
                    trace   = object$trace){
  
  tstart <- proc.time()
  if (trace>=2)
    cat("..Fitting model with scoring\n")

  f       <- object$dataRep$n - 1
  S       <- object$dataRep$S

  vccTerms <- object$intRep$vccI
  eccTerms <- object$intRep$eccI
  nrK  <- nrow(S)
  type <- object$type
  
  logL0   <- ellK(K0, S, f)

  type    <- object$type
  
  fscale  <- f #sqrt(f)
    
  theta0      <- K2theta(object, K0, scale='free')
  curr.logL   <- logL0
  curr.theta  <- theta0
  
  curr.K      <- K0
  logLeps     <- control$logLeps * abs(logL0); ##cat("Scoring: logLeps:",logLeps,"\n")
  logL.vec    <- rep(NA, maxit)
  itcount     <- 0
  
  ## Iterate here...
  ##
  if (trace>=3) cat("...Scoring start:", "logL:", curr.logL, "\n")

  repeat {
    itcount <- itcount + 1
    xxx       <- getScore(object, curr.K, scale='free')

    Sc      <- xxx$score
    DISC    <- try(solve.default(xxx$J + (Sc%*%t(Sc)/fscale), Sc))
    if (class(DISC)=='try-error'){
      cat("Error in Fisher scoring, please report...\n")
    }

    Good          <- TRUE
    stepsize      <- 1
    stephalfcount <- 0

    repeat{
      new.theta   <- curr.theta + DISC*stepsize
      new.K       <- theta2K(object, new.theta, scale='free')
      #new.K      <- theta2K2(new.theta, vccTerms, eccTerms, nrK, type, scale='free')
      #print(new.K2-new.K)
      new.logL    <- ellK(new.K, S, f)
      diff.logL   <- new.logL - curr.logL

      if (is.na(new.logL) || diff.logL<0){
        if (stephalfcount<10){
          stephalfcount <- stephalfcount + 1
          stepsize      <- stepsize / 2
          if(trace>=3) cat ("...logL:", new.logL, "outer it:",itcount,
               "New stepsize:", stepsize, "\n")
        } else {
          Good <- FALSE
          break()
        }
      } else {
        break()
      }
    }
   
    if (trace>=3) cat("...scoring iteration:", itcount, "logL:", new.logL, "dlogL",diff.logL,"\n")

    if (Good){
      curr.logL   <- new.logL
      curr.theta  <- new.theta
      curr.K      <- new.K
      logL.vec[itcount] <- curr.logL
    }

    
    if (!Good){
      break()
    } else {
      if (diff.logL < logLeps | itcount >= maxit){
        break()
      }
    }
  }
  
  if (trace>=3) cat("...Scoring iterations:", itcount, "\n")
  logL.vec <- logL.vec[!is.na(logL.vec)]

  vn <- unlistPrim(lapply(getcc(object),names))
  names(new.theta) <- vn

  ## Back to original scale
  l               <- length(getSlot(object, "vcc"))
  new.theta[1:l]  <- exp(new.theta[1:l])  
  dimnames(new.K) <- dimnames(S)
  J       <- getScore(object, curr.K, scale='original')$J

  dimnames(J) <- list(vn, vn)

  ans      <- list(K=new.K, logL=curr.logL, coef=new.theta, J=J, logL.vec=logL.vec)
  ##cat("scoring time:\n"); print(proc.time()-tstart)
  
  return(ans)
}


