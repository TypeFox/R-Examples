asNumericChar <- function(x){
##
## 1.  Convert factors to character
## 
#  print(x)
  if(length(x)<1)return(x)
  if(all(is.na(x)))return(x)
  X <- x
#  print('local copy made')
  if(is.factor(x))x <- as.character(X)
#  print('if(is.factor(x))...')
##
## 1.  Delete leading blanks and $ 
##
  x[!is.na(x)] <- tis::stripBlanks(x[!is.na(x)])
#  print(('tis::stripBlanks(x)'))
  dol <- grep('^\\$', x)
#  cat(length(dol), ' $ found: ', 
#      paste(dol, collapse=', '), '\n')
  x[dol] <- sub('^\\$', '', x[dol])
##
## 2.  find percent
##  
  pct <- grep('%$', x)
  x0 <- sub('%$', '', x)
##
## 3.  Delete commas (thousand separators) and footnote references
##
  x1 <- gsub(',', '', x0)
  x2 <- strsplit(x1, ' ')
  x. <- sapply(x2, '[', 1)
# set any blanks to NA so they don't convert to 0  
  xi <- which((!is.na(x1)) & x1=='')
#  cat(length(xi), ' blanks found: ', 
#      paste(xi, collapse=', '), '\n' )
  x.[xi] <- NA
  xo <- as.numeric(x.)
##
## 4.  rescale percents 
##
#  cat(length(pct), ' % found: ', 
#      paste(pct, collapse=', '), '\n')
  xo[pct] <- xo[pct]/100
  xo
}

asNumericDF <- function(x, keep=function(x)any(!is.na(x)),
        orderBy=NA, ignore=NULL, factors=NULL, Dates=NULL, 
        POSIX=NULL, format){
##
## 1.  Copy x
##  
  X <- as.data.frame(x)
##  
## 2.  Confirm that ignore, factors, Dates, and POSIX
##     all refer to columns of x and do not overlap.  
##
  k <- ncol(x)
  Names <- colnames(x)
#   check for Names in referenceList    
  if(is.numeric(ignore)){
    if(any(ignore<1)){
      stop('numeric ignore < 1')
    }
    if(any(ignore>k)){
      stop('ignore numeric > ncol(x)')
    }
    ignore <- colnames(x)[ignore]
  } else {
    if(length(igoops <- which(!(ignore %in% Names)))>0){
      stop('ignore = ', ignore[igoops[1]], 
           ' not in names(x) = ', 
           paste(Names, collapse=', ') )
    }
  }
# skip tests of factors, Dates, and POSIX
# ; implement later
##  
## 3.  Convert factors, Dates, and POSIX 
##
  for(f in factors){
    X[, f] <- factor(x[, f])
  }
  for(d in Dates){
    dd <- try(as.Date(x[, d], format))
    if(is(dd, 'try-error')){
      dd1 <- try(as.Date(x[, d], '%m-%d-%Y'))
      dd2 <- try(as.Date(x[, d], '%m/%d/%Y'))
      if(is(dd1, 'try-error')){
        if(is(dd2, 'try-error')){
          msg <- paste0('Failed to convert date ', 
            d, ' = ', x[1, d], ', ...')
          stop(msg)
        } else {
          X[, d] <- dd2
        } 
      } else {
        if(is(dd2, 'try-error')){
          X[, d] <- dd1 
        } else {
          na1 <- sum(is.na(dd1))
          na2 <- sum(is.na(dd2))
          if(na1<na2){
            X[, d] <- dd1
          } else {
            if(na1>na2){
              X[, d] <- dd2   
            } else {
              d1. <- sum(abs(dd1 - as.Date1970(0)), na.rm=TRUE)
              d2. <- sum(abs(dd2 - as.Date1970(0)), na.rm=TRUE)
              if(d1.<d2.){
                X[,d] <- dd1
              } else X[, d] <- dd2
            }
          }
        }
      }
    } else {
      dd1 <- try(as.Date(x[, d], '%m-%d-%Y'))
      dd2 <- try(as.Date(x[, d], '%m/%d/%Y'))
      dl <- list(dd, dd1, dd2)
      nad <- sapply(dl, function(x)sum(is.na(x)))
      naMin <- which(nad==min(nad))
      if(length(naMin)<2){
        X[, d] <- dl[[naMin]]
      } else {
        dl. <- dl[naMin]
        del <- sapply(dl., function(x){
          sum(abs(x-as.Date1970(0)), na.rm=TRUE)
        })
        delMin <- which(del==min(del))
        if(length(delMin)<1){
          X[, d] <- NA 
        } else X[, d] <- dl.[[delMin[1]]]
      }
    }
  }
  for(p in POSIX){
    if(missing(format)){
      pp <- try(as.POSIXct(x[, p]))
      if(is(pp, 'try-error')){
        msgP <- paste0('Failed to convert POSIX ', 
                      d, ' = ', x[1, d], ', ...')
        stop(msgP)
      } else X[, p] <- pp 
    } else {
      pp <- try(as.POSIXct(x[, p], format))
      if(is(pp, 'try-error')){
        msgP <- paste0('Failed to convert POSIX ', 
                       d, ' = ', x[1, d], ', ...')
        stop(msgP)
      } else X[, p] <- pp 
    }
  }
##
## 4.  Apply asNumericChar to all columns 
##     not in ignore, factors, Dates, or POSIX.  
##
  dontConvert <- union(ignore, union(factors, 
                            union(Dates, POSIX)))
  notNum <- (Names %in% dontConvert)
  numCols <- Names[!notNum]
  for(n in numCols){
    w0 <- options(warn=-1)
#    cat(colnames(x)[n], ":")
#    print(x[,n])
    xn <- asNumericChar(x[, n])
    options(warn=w0$warn)
    xnNewNA <- which(is.na(xn) & !is.na(x[, n]))
    if(length(xnNewNA)>0){
      msg0 <- paste0('NAs introduced by coercion ', 
         'in asNumericChar(c(' )
      if(length(xnNewNA)>4){
        msg1 <- paste0(msg0, paste(xnNewNA[1:4], collapse=', '), 
                       ', ...')
      } else {
        msg1 <- paste0(msg0, paste(xnNewNA, collapse=', '))
      }
      msg <- paste0(msg1, '), ', n, ')')
      warning(msg)
    }
    X[, n] <- xn
  }
##
## 5.  Keep columns specified by keep.  
##
  kp <- rep(FALSE, k)
  names(kp) <- Names
  kp[notNum] <- TRUE
  if(is.function(keep)){
      Keep <- sapply(X[numCols], keep)
      kp[numCols[Keep]] <- TRUE
  } else {
    if(!is.null(keep)) {
      if(is.logical(keep)) {
        kp[keep] <- TRUE 
      } else {
        kp[keep] <- TRUE 
      }
    }
  }
#
#  if(missing(orderBy)){
#      orderBy <- 1:length(X)
#  }
  if((length(orderBy)>1) && !is.na(orderBy)){
    o <- do.call(order, X[orderBy])
    return(X[o, kp])
  }
  X[, kp]
}

