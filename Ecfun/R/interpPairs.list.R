interpPairs.list <- function(object, .proportion, envir=list(), 
        pairs=c('1'='\\.0$', '2'='\\.1$', replace0='', 
                replace1='.2', replace2='.3'),     
        validProportion=0:1, message0=character(0), ...){
##
## 1.  Setup:  ALL.OUT, etc.    
##
# 1.1.  message0 
  if(sum(nchar(message0))<1){
    message0 <- deparse(substitute(object), 25)[1]
  }
  if((ns <- length(message0))>1){
    warning('length(message0) = ', ns, "; message0[1:2] = ", 
            paste(message0[1:2], collapse='; '), 
            ';  using the first non-null')
    message0 <- message0[nchar(message0)>0][1]
  }  
# 1.2.  check .proportion     
  if(!is.numeric(.proportion)){
    stop(message0, ': .proportion is not numeric; ', 
         ' class = ', class(.proportion))
  }  
# 1.3.  length(.proportion)>0 
  if(length(.proportion)<1){
    stop(message0, ': length(.proportion)==0')
  }
# 1.4.  is.na? 
  oops.na <- which(is.na(.proportion))
  if(length(oops.na)>0){
    stop(message0, ':  .proportion[', oops.na[1], 
         '] = NA')
  }
# 1.5.  class(validProportion)
  if(!is.numeric(validProportion)){
    stop(message0, ':  validProportion is not numeric;', 
         '  class = ', class(validProportion))
  }
# 1.6.  length(validProportion)
  if(length(validProportion)!=2){
    stop(message0, ': length(validProportion) must be 2;', 
         '  is ', length(validProportion))
  }
# 1.7.  is.na(validProportion)
  oops <- which(is.na(validProportion))
  if(length(oops)>0){
    stop(message0, ': is.na(validProportion[', oops[1], 
         '])')
  }
# 1.8.  Rows to keep 
  In <- ((validProportion[1] <= .proportion) & 
           (.proportion <= validProportion[2]) )
# 1.9.  if none, return a no op   
  if(!any(In)){
    return(list(fun='return', value=NULL))
  }  
##
## 2.  find pairs
##
  Names <- checkNames(object, 
              avoid=pairs[c(1, 4, 2, 5)])
  Chgd <- rep(FALSE, length(Names))
  names(Chgd) <- Names
#
  suf1 <- grep(pairs[1], Names, value=TRUE)
  suf2 <- grep(pairs[2], Names, value=TRUE)
  len.suf1 <- rep(NA, length(suf1))
  len.suf2 <- rep(NA, length(suf2))
  names(len.suf1) <- Names[suf1]
  names(len.suf2) <- Names[suf2]
##
## 3.  Convert identified names to common names
##
  suf1. <- sub(pairs[1], pairs[3], suf1)
  suf2. <- sub(pairs[2], pairs[3], suf2)
  suf. <- unique(c(suf1., suf2.))
  len.suf <- rep(NA, length(suf.))
  names(len.suf) <- suf. 
##
## 4.  Envir[[..]] <- eval(object) 
##
  Envir <- envir 
#  interpObj <- object
  nel <- length(object)
  Envir$.proportion <- .proportion 
#  ne <- length(Envir)
  for(j in seq(length=nel)){
#   eval(interpObj[[j]])    
#   Is this a symbol to which to assign?      
    if((j==2) && (as.character(object[[1]])=='<-')){
      Envir[[Names[j]]] <- object[[j]]
      next 
    }
    Nj <- eval(object[[j]], Envir) 
    if(is.null(Nj)){
      oj <- object[[j]]
      ojc <- deparse(oj, width.cutoff=25) 
      stop('NULL returned from eval(', Names[j], 
           ' = ', ojc, ', ...)')
    }    
    if(length(Nj)<1){
      oj <- object[[j]]
      ojc <- deparse(oj, width.cutoff=25)
      stop('eval(', Names[j], ' = ', ojc, ')\n',
           '   returned ', class(Nj), "(0)")      
    }
    Envir[[Names[j]]] <- Nj 
#   Is this a pair name? 
    s1 <- which(suf1 == Names[j])
    s2 <- which(suf2 == Names[j])
    k1 <- length(s1)
    k2 <- length(s2)
    k12 <- k1+k2
#   k12 = 0 or 1; can't be 2 
    if(k12>0){
      if(k1>0){
        len.suf1[s1] <- NROW(Nj)
#       suf1[s1];  look for match in suf2 
        j2 <- which(suf2. == suf1.[s1])
        if(length(j2)>0){
          if(suf2[j2] %in% Names[1:j]){
#         Both suf1[s1] and suf2[j2] have been eval'ed
#         Add the interpolation 
#            N12 <- (interpObj[[suf1[s1]]]*() 
#                    + interpObj[[suf2[j2]]]*proportion ) 
            argNms <- c(suf1[s1], suf2[j2], message0)
            N12 <- interpChar(Envir[c(suf1[s1], suf2[j2])], 
                      .proportion=.proportion, argnames=argNms)
            Envir[[suf2.[j2]]] <- N12 
            len.suf[suf2.[j2]] <- NROW(N12)
          } 
#         match found but not processed yet
          next
        } else {
#         match not found:
#          if(is.numeric(interpObj[[j]])){ 
#           If numeric, store interpObj[[j]] as the match 
#           with a warning     
          argNms <- c(Names[j], '.proportion', message0)
          N1 <- interpChar(x=Envir[[Names[j]]], 
                  .proportion=.proportion, argnames=argNms) 
          Envir[[suf1.[s1]]] <- N1
          len.suf[suf1.[s1]] <- NROW(N1)
        }
      } else {
        len.suf2[s2] <- NROW(Nj)
#       k2=1, because k1=0         
#       suf2[s2];  look for match in suf1 
        j1 <- which(suf1. == suf2.[s2])
        if(length(j1)>0){
          if(suf1[j1] %in% Names[1:j]){
#         Both suf1[s1] and suf2[j2] have been eval'ed
#         Add the interpolation 
#            N21 <- (interpObj[[suf1[j1]]]*(1-proportion) 
#                    + interpObj[[suf2[s2]]]*proportion ) 
            argNms <- c(suf1[j1], suf2[s2], message0)
            N21 <- interpChar(Envir[c(suf1[j1], suf2[s2])], 
                              .proportion=.proportion, argnames=argNms)
            Envir[[suf1.[j1]]] <- N21 
            len.suf[suf1.[j1]] <- NROW(N21)
          } 
#         match found but not processed yet
          next
        } else {
#         match not found;  store interpObj[[j]] as the match           
          argNms <- c(Names[j], '.proportion', message0) 
          N2 <- interpChar(Nj, 
                        .proportion=.proportion, argnames=argNms) 
          Envir[[suf2.[s2]]] <- N2
          len.suf[suf2.[s2]] <-NROW(N2)
          next 
        }     
      }
    }
  }
##
## 5.  Other vectors or data.frames 
##     with the same number of rows?  
##
  Drop <- (Names %in% c(suf1, suf2))
  Keep. <- c(Names[!Drop], suf.)
  interpOut <- Envir[Keep.]  
  if(length(interpOut)>0){
    objLen <- sapply(interpOut, NROW)
  } else objLen <- integer(0)
#
  nSuf <- length(suf.)
  if(nSuf>0){
    ln <- max(objLen[suf.])
  } else ln <- (-Inf)
#  lp <- length(.proportion)
  lp <- length(In)
  if(lp < ln) {
#    .proportion <- rep(.proportion, length=ln)
    In <- rep(In, length=ln)
    N <- ln
  } else {      
    N <- lp 
    if((lp>ln) && (ln>1)){
      msg <- paste0('length(.proportion) = ', lp, 
              ' > max length(pairs) = ', ln, 
              ';\n  pairs found for ', 
              paste(suf., collapse=', '))
      warning(msg)
    }
  }
# Cols to trim? 
  cols2trim <- ((objLen==N) & (N>1))
# trim   
  for(s. in names(interpOut)[cols2trim]){
#   Retain only "In" in s.
    S. <- interpOut[[s.]]
    sub. <- any(is.language(S.) | is.function(S.))
    if(sub.)next 
#    
    ndim <- length(dim(S.))
    if(ndim<2){
      Subset <- ((!is.null(S.)) && !is.function(S.))
      if(Subset){
          S. <- S.[In]
      }
    } else {     
      sub2 <- (is.data.frame(S.) || (is.matrix(S.)))
      if(sub2){
          S. <- S.[In,, drop=FALSE]
      } 
    }
    interpOut[[s.]] <- S.
    Chgd[s.] <- TRUE 
  }
##
## 6.  Restore object[!Chgd]
## 
  Nam3 <- names(Chgd)
  Nout <- (Nam3 %in% names(interpOut))
  noChg <- (Nout & !Chgd)  
  noC <- Nam3[noChg]
  for(nc in noC){
#   in case is.null(names(object))
    nc. <- which(Names==nc)
    interpOut[[nc]] <- object[[nc.]]
  }
  interpOut
}