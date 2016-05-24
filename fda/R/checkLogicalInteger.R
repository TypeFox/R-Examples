checkLogical <- function(x, length., warnOnly=FALSE) {
##
## 0.  Set up
##
  onExit <- function(...){
    Msg <- paste("In ", callEnv, ':  ', ..., sep='') 
    if(warnOnly){
      warning(Msg, call.=FALSE)
      return(FALSE)
    }
    else
      stop(Msg, call.=FALSE)
  }
##
## 1.  External name of 'x'
##  
  xName <- substring(deparse(substitute(x)), 1, 44)
# name of calling function
  callEnv <- sys.call(-1)[[1]]
  if(is.null(callEnv)) callEnv <- sys.call()[[1]]
##
## 2.  is.null(x)   
##
  good <- TRUE
  if(is.null(x))
    good <- onExit('is.null(', xName, ')')
# good <- stopWarn('is.null(', xName, ')', warnOnly=warnOnly, which.parent=which.parent+1)
##
## 3.  check class(x)
##  
  if(class(x) != 'logical')
    good <- (good & onExit('class(', xName, ') = ', class(x),
                   ";  should be 'logical'")) 
##
## 4.  Check 'length.' 
##
  if(!missing(length.)){
    if(length(x) != length.){
      good <- (good & onExit('length(', xName, ') = ', length(x), ' != ',
             ' required length = ', length.) ) 
    }
  }
##
## 5.  Done
##  
  good
}

checkNumeric <- function(x, lower, upper, length., integer=TRUE,
           unique=TRUE, inclusion=c(TRUE,TRUE), warnOnly=FALSE){
##
## 0.  set up exit / return processing with warnOnly
##
  onExit <- function(...){
    Msg <- paste("In ", callEnv, ':  ', ..., sep='') 
    if(warnOnly){
      warning(Msg, call.=FALSE)
      return(FALSE)
    }
    else
      stop(Msg, call.=FALSE)
  }
##
## 1.  External name of 'x'
##  
  xName <- substring(deparse(substitute(x)), 1, 44)
# name of calling function
  callEnv <- sys.call(-1)[[1]]
  if(is.null(callEnv)) callEnv <- sys.call()[[1]]
##
## 2.  is.null(x)   
##
  if(is.null(x))return(TRUE)  
##
## 3.  check class(x)
##
  if(!is.numeric(x))
    onExit('class(', xName, ') = ', class(x), ";  should be 'numeric'")
  if(integer){
    x. <- round(x)
    d <- abs(x.-x)
    d.ne0 <- (d != 0)
    if(any(d.ne0)){
      id <- which(d == max(d))[1]
      nne <- sum(d.ne0)
      if(nne>1)
        onExit(sum(d.ne0), ' non-integer values;  ',
               'the most extreme is ', xName, '[', id, '] = ', x[id])
      else
        onExit('One non-integer value:  ',
               xName, '[', id, '] = ', x[id])
    }
  }
##
## 4.  Check limits 
##
  dLow <- (x-lower)
  xLow <- {
    if(inclusion[1]) (dLow<0) else (dLow<=0) 
  }
  if(any(xLow)){
    ilo <- which(dLow == min(dLow))[1]
    nlo <- sum(xLow)
    if(nlo>1)
      onExit(nlo, ' low values;  the most extreme is ',
             xName, '[', ilo, '] = ', x[ilo])
    else
      onExit('One low value:  ', xName, '[', ilo, '] = ', x[ilo])
  }
#
  dHi <- (x-upper)
  inclusion <- rep(inclusion, length=2)
  xHi <- {
    if(inclusion[2]) (dHi>0) else (dHi>=0)
  }
  if(any(xHi)){
    ihi <- which(dHi==max(dHi))[1]
    nhi <- sum(xHi)
    if(nhi>1) 
      onExit(nhi, ' high values;  the most extreme is ',
                   xName, '[', ihi, '] = ', x[ihi])
    else
      onExit('One high value:  ', xName, '[', ihi,
                   '] = ', x[ihi]) 
  }
##
## 5.  Check unique
##
  if(length(x)>1){
    x. <- sort(x)
    dx <- diff(x.)
    if(any(dx==0)){
      x <- unique(x) 
      iun <- which(dx==0)[1]
      nun <- sum(dx==0)
      if(nun>1) 
        onExit(nun, ' repeated values in ', xName,
                     ';  the smallest is ', x.[iun])
      else
        onExit('One repeated value in ', xName, ':  ', x.[iun])
    }
  }
##
## 6.  Check length 
##
  if(!missing(length.)){
    if(length(x) != length.)
      onExit('length(', xName, ') = ', length(x), ' != ',
             ' required length = ', length.)
  }
##
## 7.  Done
##    
  x 
}

checkLogicalInteger <- function(x, length., warnOnly=FALSE){
##
## 0.  set up exit / return processing with warnOnly
##
  onExit <- function(...){
    Msg <- paste("In ", callEnv, ':  ', ..., sep='') 
    if(warnOnly){
      warning(Msg, call.=FALSE)
      return(FALSE)
    }
    else
      stop(Msg, call.=FALSE)
  }
##
## 1.  External name of 'x'
##  
  xName <- substring(deparse(substitute(x)), 1, 44)
# name of calling function
  callEnv <- sys.call(-1)[[1]]
  if(is.null(callEnv)) callEnv <- sys.call()[[1]]
##
## 2.  is.null(x)   
##
  if(is.null(x))return(rep(TRUE, length=length.))  
##
## 3.  check class(x)
##
#  3.1.  is.logical?    
  if(is.logical(x)){
    if(missing(length.)) return(x)
    {
      if(length(x) == length.) return(x)
      else
        onExit('length(x) = ', length(x), ';  should be ',length.)
    }
  }
#  3.2.  is.numeric?  
  if(!is.numeric(x))
    onExit('class(', xName, ') = ', class(x),
           ";  should be 'numeric' or 'logical'")
#  3.3.  is.integer?  
  x. <- round(x)
  d <- abs(x.-x)
  d.ne0 <- (d != 0) 
  if(any(d.ne0)){
    id <- which(d == max(d))[1]
    nne <- sum(d.ne0)
    if(nne>1)
      onExit(sum(d.ne0), ' non-integer values;  ',
             'the most extreme is ', xName, '[', id, '] = ', x[id])
    else
      onExit('One non-integer value:  ',
             xName, '[', id, '] = ', x[id])
  }
##
## 4.  Check limits 
##
#  4.1.  x < 1?    
  xLow <- (x<1) 
  if(any(xLow)){
    ilo <- which(x == min(x))[1]
    nlo <- sum(xLow)
    if(nlo>1)
      onExit(nlo, ' low values;  the most extreme is ',
             xName, '[', ilo, '] = ', x[ilo])
    else
      onExit('One low value:  ', xName, '[', ilo, '] = ', x[ilo])
  }
#   4.2.  x > length.?
  if(missing(length.) && warnOnly){
    onExit("argument 'length.' is missing;  setting to max(x)")
    length. <- max(x) 
  }
  xHi <- (x>length.) 
  if(any(xHi)){
    ihi <- which(x==max(x))[1]
    nhi <- sum(xHi)
    if(nhi>1) 
      onExit(nhi, ' high values;  the most extreme is ',
                   xName, '[', ihi, '] = ', x[ihi])
    else
      onExit('One high value:  ', xName, '[', ihi,
                   '] = ', x[ihi]) 
  }
##
## 5.  Check unique
##
  if(length(x)>1){
    x. <- sort(x)
    dx <- diff(x.)
    if(any(dx==0)){
      x <- x[dx != 0] 
      iun <- which(dx==0)[1]
      nun <- sum(dx==0)
      if(nun>1) 
        onExit(nun, ' repeated values in ', xName,
                     ';  the smallest is ', x.[iun])
      else
        onExit('One repeated value in ', xName, ':  ', x.[iun])
    }
  }
##
## 6.  Convert to logical  
##
  X <- rep(FALSE, length.)
  X[x] <- TRUE 
##
## 7.  Done
##    
  X
}
  

