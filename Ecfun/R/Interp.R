Interp <- function(x, ...){
  UseMethod('Interp')
}

InterpChkArgs <- function(x, y, proportion, 
        argnames=character(3), message0=character(0), ...){
##
## 1.  proportion missing or not numeric
##
  nmx <- 'x'
  nmy <- 'y'
  nmp <- 'proportion'
  if(nchar(argnames[1])>0){
    nmx <- paste0(nmx, ' (', argnames[1], ')')
  }
  if(nchar(argnames[2])>0){
    nmy <- paste0(nmy, ' (', argnames[2], ')')
  }
  if(nchar(argnames[3])>0){
    nmp <- paste0(nmp, ' (', argnames[3], ')')
  }
#  
  if(missing(proportion)){
    msg0.p <- paste0(message0, ":  ", nmp, ' is missing.')
    stop(msg0.p)
  }
  if(!is.numeric(proportion)){
    msg1.p <- paste0(message0, ":  ", nmp, 
                     ' is not numeric.')
    stop(msg1.p)
  }
##
## 2.  missing x 
##
  np <- length(proportion)
  if(missing0(x)){ 
    if(missing0(y)){
      msg0.xy <- paste0(message, ":  ", nmx, 
              " and ", nmy, ' are both missing.')
      stop(msg0.xy)      
    }
    compyp <- compareLengths(y, proportion, name.x=argnames[2], 
                        name.y=argnames[3], message0=message0)
    ny <- length(y)
    outList <- list(xout=rep(y, length=max(np, ny))) 
  } else {
##
## 3.  x provided.  y missing?   
##
    if(missing0(y)){
      compxp <- compareLengths(x, proportion, name.x=argnames[1], 
                               name.y=argnames[3], message0=message0)
      nx <- length(x)
      outList <- list(xout=rep(x, length=max(nx, np)))
    } else {
##
## 4.  both x and y provided.  
##      
      nx <- length(x)
      ny <- length(y)
      if(all(proportion<=0)){
        compxp <- compareLengths(x, proportion, name.x=argnames[1], 
                                 name.y=argnames[3], message0=message0)
        outList <- list(xout=rep(x, length=max(nx, np)))
      } else if(all(proportion>=1)){
        compyp <- compareLengths(y, proportion, name.x=argnames[2], 
                                 name.y=argnames[3], message0=message0)
        outList <- list(xout=rep(y, length=max(np, ny))) 
      } else {
        compxy <- compareLengths(x, y, name.x=argnames[1], 
            name.y=argnames[2], message0=message0)
        compxp <- compareLengths(x, proportion, name.x=argnames[1], 
            name.y=argnames[3], message0=message0)
        compyp <- compareLengths(y, proportion, name.x=argnames[2], 
            name.y=argnames[3], message0=message0)
        #
        pLength1 <- (np==1)
##
## 5.  Character or Numeric?
##
# Character if: 
# (1) At least one of x and y is character.  
# (2) At least one of x and y is neither 
#     logical, integer, numeric, complex nor raw, 
#     and class(unclass(.)) is either integer or character.  
        clx <- class(x)
        cly <- class(y)
        cix <- classIndex(x)
        ciy <- classIndex(y)
        ux <- unclass(x)
        uy <- unclass(y)
        clux <- class(ux)
        cluy <- class(uy)
        #        
#        basicNum <- c('logical', 'integer', 'numeric', 
#                      'complex', 'raw')
#        basicNumx <- all(clx %in% basicNum) 
#        basicNumy <- all(cly %in% basicNum)
        basicNumx <- (cix < 7)
        basicNumy <- (ciy < 7)
        #
        unclassCh <- c('integer', 'character')
        unclassChx <- (clux %in% unclassCh)
        unclassChy <- (cluy %in% unclassCh)
        #
        Ch2x <- ((!basicNumx) & unclassChx)
        Ch2y <- ((!basicNumy) & unclassChx)
        #
        Ch <- (((!basicNumx) | (!basicNumy)) & 
               (Ch2x | Ch2y) )
        # 
        Raw <- ((!Ch) & is.raw(clux) & is.raw(cluy))
        #
        if(!basicNumy){
          outClass <- attributes(y)
        } else {
          if(!basicNumx){
            outClass <- attributes(x)
          } else outClass <- NA 
        }
# Otherwise, 
# the list includes components "algorithm", "x", "y", 
# "proportion", "pLength1" (defined below), "raw", and 
# "outclass".  The "algorithm" component must be either 
# "Numeric" or "Character".
        nxyp <- max(nx, ny, np)
        P <- rep(proportion, length=nxyp)
        if(Ch){
          X <- rep(as.character(ux), length=nxyp)
          Y <- rep(as.character(uy), length=nxyp)
          outList <- list(algorithm='Character', 
              x=X, y=Y, proportion=P, pLength1=pLength1, 
              raw=Raw, outClass=outClass)                    
        } else {
          X <- rep(as.numeric(ux), length=nxyp)
          Y <- rep(as.numeric(uy), length=nxyp)
          outList <- list(algorithm='Numeric', 
              x=X, y=Y, proportion=P, pLength1=pLength1, 
              raw=Raw, outClass=outClass)                              
        }        
      }
      
    }
  }
  outList 
}  
    
Interp.default <- function(x, y, proportion, 
            argnames=character(3), message0=character(0), ...){
##
## 1.  InterpChkArgs 
##  
  argsChk <- InterpChkArgs(x, y, proportion, 
            argnames, message0, ...)
##
## 2.  Done?  
##
  if('xout' %in% names(argsChk)){
    xout <- argsChk$xout
  } else if(argsChk$algorithm=='Numeric'){
    xout <- InterpNum(argsChk, ...)
  } else if(argsChk$algorithm=='Character'){
    xout <- InterpChar(argsChk, ...)
  } else {
    print(argsChk)
    print(message0)
    print(argnames)
    stop('bug(?) in InterpChkArgs()$algorithm = ', 
         argsChk$algorithm, ' is neither ', 
         ' Numeric nor Character')
  }
  xout
}

InterpNum <- function(argsChk, ...){
  xout <- with(argsChk, (1-proportion)*x + proportion*y)
  if(argsChk$raw){
    xout <- as.raw(xout)
  }
  if(!is.na(argsChk$outClass)){
    oC <- argsChk$outClass
    oC$.Data <- xout
    xout <- do.call(structure, oC)
  }
  xout  
}

InterpChar <- function(argsChk, ...){
## 
## 1.  Base on x or y?  
## 
  nch.x <- nchar(argsChk$x)
  nch.y <- nchar(argsChk$y) 
  P <- argsChk$proportion 
  N <- length(argsChk$x)
#
  swap <- (nch.y<nch.x)
  Ny <- nch.y
  Nx <- nch.x
  Z <- argsChk$y
  Z[swap] <- argsChk$x[swap]
  Ny[swap] <- nch.x[swap]
  Nx[swap] <- nch.y[swap]
  P[swap] <- (1-P[swap])
#
  dxy <- (Ny-Nx)
##
## 2.  pLength1 
##
  if(argsChk$pLength1){
    Dxy <- cumsum(dxy)
    Nch.Dxy <- tail(Dxy, 1)
    Nch1 <- (round(P*Nch.Dxy)-c(0, head(Dxy, -1)))
  } else {
    Nch1 <- (Nx+round(P*dxy))
  }
##
## 3.  substring 
##
#  keepCh <- (nch.x+pmax(0, Nch1-Dxy))    
#  Out <- substring(Z, 1, keepCh)
  Out <- substring(Z, 1, Nch1)
##
## 4.  Done 
## 
  Out 
}     
