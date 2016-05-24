interpChar <- function(x, ...){
  UseMethod('interpChar')
}

interpChar.list <- function(x, .proportion, 
         argnames=character(3), message0=character(0),  ...){
##
## 1.  message0?  
##
  if(sum(nchar(message0))==0){
    message0 <- deparse(substitute(x), 25)
  }
##
## 2.  length(x)<2
##
  if(length(x)<2){
#    lx <- length(x[[1]])
#    lp <- length(.proportion)
    xNm <- names(x)
    name0 <- FALSE 
    name.x <- argnames[1]
    if(nchar0(x)) {
      name0 <- TRUE
      name.x <- createMessage(xNm, 35, default='x') 
    }
    if(nchar0(xNm)) xNm <- name.x 
    name.y <- argnames[2] 
    if(nchar0(name.y)) {
      name.y <- '.proportion'
      name0 <- TRUE 
    }
    message0 <- paste(message0, argnames[3]) 
    if(name0){
      message0 <- paste0('in interpChar.list:', message0)
    }  
    compareLengths(x[[1]], .proportion, name.x, name.y,
                   message0, ...)    
    if(is.numeric(x[[1]])){
      if(nchar0(xNm)){
        warning('numerical interpolation in a list of length 1', 
                '\n returns the input')
      } else {
        warning('numerical interpolation in a list of length 1', 
                '\n with an element named ', xNm, 
                ';  returning the input.')
      }
      return(x[[1]])
    }
    out <- interpChar.default('', x[[1]], .proportion, ...)
    return(out)
  }
##
## 3.  length(x)>1 
##
  interpChar.default(x[[1]], x[[2]], .proportion, 
                     argnames, message0, ...)
}

interpChar.default <- function(x, y, .proportion, 
           argnames=character(3), message0=character(0), ...){
##
## 1.  message0?  
##
  misx <- missing(x)
  misy <- missing(y)
  if(sum(nchar(message0))==0){
    message0 <- createMessage(deparse(x), 25L)
  }
##
## 2.  compareLengths(x, .proportion, ...)  
##  
  message0 <- paste(message0, argnames[3])
  name.x <- argnames[1]
  name0 <- FALSE 
  if(nchar0(name.x)) {
    message0 <- paste0('in interpChar.default:', message0)
    name.x <- 'x'
    name0 <- TRUE 
  }
  name.p <- '.proportion'
#
  if(misx || is.null(x)){ 
    if(misy || is.null(y)){
      warning(message0, ':  both x and y are missing or NULL;', 
              '   returning NULL')
      return(NULL)
    }
    ciy <- classIndex(y)
    if(misx){ 
      x <- do.call(index2class(ciy), list(length(y)))
    } else {
      x <- createX2matchY(x, y)
    }
  } else {
    if(misy || is.null(y)){
      cix <- classIndex(x)
      if(misy){
        y <- do.call(index2class(cix), list(length(x)))      
      } else {
        y <- createX2matchY(y, x)
      }
    }
  }
  nx <- length(x)
  ny <- length(y)
  cix <- classIndex(x)
  ciy <- classIndex(y)
  if(nx<1){
    if(ny<1){
      return(createX2matchY(x, y))
    }
    x <- createX2matchY(x, y)    
    if(cix>ciy)
      y <- as(y, index2class(cix))
  } else {
    if(ny<1)
      y <- createX2matchY(y, x)
  }
  cL <- compareLengths(x, .proportion, name.x, name.p,
                 message0, ...)    
##
## 3.  numeric? 
##  
#  if(missing(y)){
#  3.1.  missing(y)    
#    if(is.numeric(x)){
#      warning('numeric interpolation with one input;', 
#              '  returning that.')
#      return(x)
#    }
#    y <- x
#    x <- '' 
#  } else { 
#  3.2.  y is not missing:  Check lengths     
  name.y <- argnames[2]
  if(is.null(name.y) || (nchar(name.y)==0)){
      name.y <- 'y'
      if(!name0){
        message0 <- paste0('in interpChar.default:', message0)
      }
  }  
  cL.xy <- compareLengths(x, y, name.x, name.y,
                   message0, ...)    
#   Numeric?      
  ciz <- max(cix, ciy)
  if(ciz<6){
      out <- (x*(1-.proportion) + y*.proportion) 
      return(out)
  }
##
## 4.  not numeric
##
#  4.1.  as.character 
  xc <- as.character(x)
  yc <- as.character(y)
#  4.2.  Same length
#  nx <- length(xc)
#  ny <- length(yc)
  np <- length(.proportion)
  N <- max(nx, ny, np)
  X <- rep(xc, length=N)
  Y <- rep(yc, length=N)
  P. <- rep(.proportion, length=N)
#  2.3.  number of characters 
  nch.y <- nchar(Y)
  nch.x <- nchar(X)
  swap <- (nch.y<nch.x)
  Z <- Y
  Z[swap] <- X[swap]
  Ny <- nch.y
  Nx <- nch.x
  P <- P. 
  Ny[swap] <- nch.x[swap]
  Nx[swap] <- nch.y[swap]
  P[swap] <- (1-P.[swap])
  dxy <- (Ny-Nx)
  Dxy <- cumsum(dxy)
  DN <- Dxy[N]
  cumCh <- P*DN
# 0 <= DN < Inf, 
# so the only way is.na(cumCh) 
# is abs(P)==Inf & DN==0 
  cumCh[is.na(cumCh)] <- P[is.na(cumCh)]
#  
  D.xy <- c(0, Dxy[-N])
  Pd <- (cumCh - D.xy) 
#  Pd[is.na(Pd)] <- 1 
#  Pd. <- pmax(0, pmin(Pd, 1))
  Out <- character(N)
  sely <- round(Nx + Pd)
  Out[sely>0] <- substring(Z[sely>0], 1, sely[sely>0])
##
## 4.  Done 
## 
  Out 
}