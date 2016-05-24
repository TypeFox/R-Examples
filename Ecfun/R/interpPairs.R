interpPairs <- function(object, ...){
  if(is.null(object)){
    warning('interpPairs:  object = NULL.  Possible error;', 
            '  returning NULL')
    return(NULL)
  }
  if(is(object, '<-')){
    return(interpPairs.function(object, ...))
  }
#  print(str(object))
#  print(deparse(substitute(object)))
  UseMethod('interpPairs')  
}

#"interpPairs.<-" <- function(object, 
#
# *** THIS DID NOT WORK:  
# *** S3 methods dispatch seemed not to recognize "<-" as a class.
# *** It does use "call";  see interpPairs.call
#
#     nFrames=1, iFrame=nFrames, endFrames=round(0.2*nFrames), 
#     envir = parent.frame(), 
#     pairs=c('1'='\\.0$', '2'='\\.1$', replace0='', 
#                               replace1='.2', replace2='.3'),     
#     validProportion=0:1, message0=character(0), ...){
  #  co <- class(object)  
  #  class(object) <- 'function' 
  #  ip <- nextMethod()???
  #  class(ip) <- co
  #  return(ip) ???  
#  interpPairs.function(object, nFrames=nFrames, 
#      iFrame=iFrame, endFrames=endFrames, envir=envir, 
#      pairs=pairs, validProportion=validProportion, 
#      message0=message0, ...)
#}

interpPairs.call <- function(object, 
            nFrames=1, iFrame=nFrames, 
            endFrames=round(0.2*nFrames), 
            envir = parent.frame(), 
            pairs=c('1'='\\.0$', '2'='\\.1$', replace0='', 
                     replace1='.2', replace2='.3'),     
            validProportion=0:1, message0=character(0), ...){
#  co <- class(object)  
#  class(object) <- 'function' 
#  ip <- nextMethod()???
#  class(ip) <- co
#  return(ip) ???  
  interpPairs.function(object, 
        nFrames=nFrames, iFrame=iFrame, endFrames=endFrames, 
        envir=envir, pairs=pairs, 
        validProportion=validProportion, message0=message0, ...)
}

interpPairs.function <- function(object, 
              nFrames=1, iFrame=nFrames, 
              endFrames=round(0.2*nFrames), 
              envir = parent.frame(), 
              pairs=c('1'='\\.0$', '2'='\\.1$', replace0='', 
                       replace1='.2', replace2='.3'),     
              validProportion=0:1, message0=character(0), ...){
##
## 1.  Setup:  ALL.OUT, etc.    
##
# 1.1.  message0 
  if(sum(nchar(message0))<1){
    message0 <- createMessage(deparse(object))
  }
  if((ns <- length(message0))>1){
    warning('length(message0) = ', ns, "; message0[1:2] = ", 
            paste(message0[1:2], collapse='; '), 
            ';  using the first non-null')
    message0 <- message0[nchar(message0)>0][1]
  }  
##
## 2.  look for firstFrame, lastFrame, and Keep
##
  Names <- names(object)
  Object <- object 
  Envir <- new.env() 
  if('firstFrame' %in% Names){
    firstF <- eval(object$firstFrame, Envir, envir) 
    if(is.null(firstF))firstF <- 1
    Object$firstFrame <- NULL 
  } else firstF <- 1
  Envir$firstFrame <- firstF 
#
  if(length(nFrames)<1){
    nFrames <- firstF 
  }
  if('lastFrame' %in% Names){
    lastF <- eval(object$lastFrame, Envir, envir)  
    if(is.null(lastF)) lastF <- (nFrames-endFrames+1) 
    Object$lastFrame <- NULL 
  } else lastF <- (nFrames-endFrames+1) 
#
  Envir$lastFrame <- lastF 
#
  if('Keep' %in% Names){
    Keep <- eval(object$Keep, Envir, envir)  
    Object$Keep <- NULL 
  } else Keep <- TRUE 
  Envir$Keep <- Keep 
##
## 3.  .proportion 
##
  chk.lf <- compareLengths(lastF, firstF, 
          name.x='lastFrame', name.y='firstFrame', 
          message0=message0)
  dF <- (lastF-firstF)
  chk.if <- compareLengths(iFrame, firstF, 
          name.x='iFrame', name.y='firstFrame', 
          message0=message0)
  pDone <- pmin((iFrame-firstF) / dF, 1)
  pDone[is.na(pDone)] <- 1
  chk.kf <- compareLengths(Keep, firstF, 
          name.x='Keep', name.y='FirstFrame', 
          message0=message0 )
  pDone[iFrame>lastF] <- (2*Keep - 1)
##
## 4.  validProportion 
##
  In <- ((validProportion[1]<=pDone) & 
           (pDone<=validProportion[2]) )
  if(!any(In)){
    return(enquote(NULL))
  }
##
## 5.  find pairs
## 
  Nam2 <- names(Object)
#  Chgd <- rep(FALSE, length(Nam2))
#  names(Chgd) <- Nam2
# suf1, suf2 = indices in Nam2 
  suf1 <- grep(pairs[1], Nam2, value=TRUE)
  suf2 <- grep(pairs[2], Nam2, value=TRUE)
  suf12 <- c(suf1, suf2)
  n1 <- length(suf1)
  n2 <- length(suf2)
  n12 <- length(suf12)
  Suf1 <- vector('list', n1)
  Suf2 <- vector('list', n2)
  Suf12 <- vector('list', n12)
  i.1 <- seq_len(n1)
# suf1., suf2. = Nam2[c(suf1, suf2)]
  suf1. <- sub(pairs[1], pairs[3], suf1)
  suf2. <- sub(pairs[2], pairs[3], suf2)
  suf. <- unique(c(suf1., suf2.))
##
## 6.  eval 
##
  inam <- seq_len(length(Nam2))
  for(i. in inam){
    si <- Nam2[i.]
#  6.1  match?      
    if(!(si %in% suf12)) next 
#  6.2.  eval si  
    Si <- eval(object[[si]], Envir, envir) 
    if(is.null(Si))Si <- logical(0)
    Envir[[si]] <- Si 
#  6.3.  match 
    k1 <- which(suf1 == si)
    if(length(k1)>0){
      s. <- suf1.[k1]
      i2 <- which(suf2. == s.)
      if(length(i2)<1){
#  6.4.  only s1 
        argnms <- c(si, '(no match)', '.proportion')
        S. <- interpChar(x=Si, y=logical(0), pDone, 
                           argnms, message0)
        if(is.null(S.))S. <- logical(0)
#        S. <- Si        
        Object[[s.]] <- S. 
        Envir[[s.]] <- S. 
        next 
      }
#  6.5. s1 and s2;  s2 yet?           
      s2 <- suf2[i2]
      if(s2 %in% Nam2[1:i.]){
        argnms <- c(si, s2, '.proportion') 
        S. <- interpChar(Si, Envir[[s2]], pDone, 
                           argnms, message0) 
        if(is.null(S.))S. <- logical(0)
        Object[[s.]] <- S. 
        Envir[[s.]] <- S. 
      } 
      next 
    }       
#  6.6.  in suf2 
    k2 <- which(suf2 == si)
#  6.7.  also suf1?  
    s. <- suf2.[k2]
    i1 <- which(suf1. == s.)
    if(length(i1)<1){
#  6.8.  only s2 
      argnms <- c('(no match)', si, '.proportion')
      S. <- interpChar(x=logical(0), y=Si, pDone, 
                         argnms, message0)
      if(is.null(S.))S. <- logical(0)
      S. <- Si 
      Object[[s.]] <- S. 
      Envir[[s.]] <- S. 
      next 
    }
#  6.9. s1 and s2;  s1 yet?           
    s1 <- suf1[i1]
    if(s1 %in% Nam2[1:i.]){
      argnms <- c(s1, si, '.proportion') 
      S. <- interpChar(Envir[[s1]], Si, pDone, 
                   argnms, message0) 
      if(is.null(S.))S. <- logical(0)
      Object[[s.]] <- S. 
      Envir[[s.]] <- S. 
    }       
  }
##
## 7.  drop 
##
  Nam3 <- names(Object)
  if(length(Nam3)==length(Object)){
    drop <- (Nam3 %in% suf12)
    Obj <- Object[!drop]
  } else Obj <- Object 
##
## 8.  common length 
##
  len.p <- length(pDone)
  if(length(suf.)>0){
    len1 <- sapply(Obj[suf.], NROW)
  } else len1 <- integer(0)
  N <- max(1, len.p, len1)
  incomp.p <- (N %% len.p)
  if(incomp.p>0){
    msg.p <- paste0(message0, ': incompatible lengths.', 
                  ' length(.proportion) = ', len.p, 
                  ';  common length = ', N) 
    warning(msg.p)
  }
  incomp1 <- which((N %% len1) > 0) 
  if(length(incomp1)>0){
    msg1 <- paste0(message0, ': incompatible lengths.', 
                 ' length(', suf.[incomp1[1]], 
                 ') = ', len1[incomp1[1]], 
                 ';  common length = ', N) 
    warning(msg1)
  }
#  pD <- rep_len(pDone, N)
  In. <- rep_len(In, N)
  for(s. in suf.){
    Obj[[s.]] <- Obj[[s.]][In.]
  }    
##
## 9.  Done 
##
  Obj 
}