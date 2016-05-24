matchName <- function(x, data, Names=1:2, 
                      nicknames=matrix(character(0), 0, 2), 
                      namesNotFound="attr.replacement", ...){
##
## 1.  Check x  
##
  x0 <- x
  if(is.data.frame(x)){
    x <- as.matrix(x)
  }
  if(!is.character(x)){
    if(is.data.frame(x)){
      stop('class(as.matrix(x)) must be character;',
           ' is ', class(x))
    }
    if(is.null(x))return(NULL)
    stop('class(x) = ', class(x), '; must be either character', 
         ' or a data.frame that can be converted to character')
  }
  dimx <- dim(x)
  if(length(dimx)<2){
    x <- parseName(x, namesNotFound=namesNotFound, ...)
  } else if(length(dimx)>2){
    stop('length(dim(x)) must be 1 or 2;  is ', length(dimx))
  }
##
## 2.  Check data
##
  dimd <- dim(data)
  lend <- length(dimd)
  if(lend != 2){
    erd <- paste0('length(dim(data)) must be 2;  is ', lend)
    if(lend>0){
      erd <- paste0(erd, "; dim(data) = ", 
                    paste(dimd, collapse=', '))
    }
    stop(erd)    
  }
##
## 3.  Names
##
  Names0 <- Names
  if(is.numeric(Names)){
    if(any(Names<1)){
      stop('Names is numeric so must be at least 1;  is ', 
           paste(Names, collapse=', '))
    }
    if(any(trunc(Names)>dimd[2])){
      stop('Names is numeric so must be at most dim(data)[2] = ', 
        dimd[2], ';  is ', paste(Names, collapse=', '))
    }
    dn <- data[, Names]
    Names <- subNonStandardNames(dn, namesNotFound=namesNotFound, 
                                 ...)
  } else if(is.logical(Names)){
    if(length(Names) != dimd[2]){
      stop('Names is logical, so its length must match dim(data)[2];', 
           '  length(Names) = ', length(Names))
    }
    dn <- data[, Names]
    Names <- subNonStandardNames(dn, namesNotFound=namesNotFound, 
                                 ...)
  } else if(is.character(Names)){
    nNms <- NROW(Names)
    if(nNms != dimd[1]){
      chk <- (Names %in% colnames(data))
      if(any(!chk)){
        stop('Names is character, not found in colnames(data): ', 
             paste(Names[!chk], collapse=', '))
      }
      dn <- data[, Names]
      Names <- subNonStandardNames(dn, namesNotFound=namesNotFound, 
                                   ...)
    }
  } else stop('class(Names) must be either numeric, logical, or ', 
              'character;  is ', class(Names))
##
## 4.  nicknames 
##
  if(!is.character(nicknames)){
    stop('nicknames must be character;  are ', class(nicknames))
  }
  dimn <- dim(nicknames)
  if(length(dimn) != 2){
    ern <- paste0('nicknames must be a matrix; length(dim(nicknames)) = ', 
                  length(dimn))
    stop(ern)
  }
##
## 5.  xlist <- matchName1(x[, 1], data, Names[, 1], data, ...)
##
  xlist <- matchName1(x[, 1], data, Names[, 1], ...)
# names(xlist) <- paste(x[, 1], x[, 2], sep=', ')
  xRestore <- paste(x[, 1], x[, 2], sep=', ')
  xR <- sub('^, ', '', xRestore)
  names(xlist) <- xR 
##
## 6.  For any component i of xlist with 0 or multiple rows 
##  
  id <- seq(length=dimd[1])
  nx <- NROW(x)
  out <- vector(nx, mode='list')
  names(out) <- names(xlist)
  for(ix in seq(length=nx)){
    jd <- xlist[[ix]]
    ni <- length(jd)
    if(ni != 1){
#      if(ni<1)jd <- id 
#   If ni < 1:  No match for first name;  next       
      if(ni<1)next 
      xi <- strsplit(x[ix, 2], ' ')[[1]] 
      for(j in seq(along=xi)){
#        jd2 <- matchName1(xi[j], data[jd,, drop=FALSE], 
#          Names[jd, -1, drop=FALSE], nicknames=nicknames, ...)[[1]]
#        jd2 <- matchName1(xi[j], data[jd,], 
#                          Names[jd, -1, drop=FALSE], nicknames=nicknames, ...)[[1]]
        jd2 <- matchName1(xi[j], data[jd,], Names[jd, -1], 
                          nicknames=nicknames, ...)[[1]]
        if(length(jd2)>0){
          jd <- jd[jd2]
#          if(length(jd)<2) next 
          if(length(jd)<2) break # j  
        } # else if(length(jd2)<2) break         
      }
    } 
    ni <- length(jd)
    if((0<ni) & (ni < dimd[1]))out[[ix]] <- data[jd,, drop=FALSE]
  }    
##
## 7.  done 
##
  notF <- c(attr(x, 'namesNotFound'), 
            attr(Names, 'namesNotFound'))
  if(length(notF)>0){
    attr(out, 'namesNotFound') <- notF
  }
  out
}