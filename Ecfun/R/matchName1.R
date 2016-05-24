matchName1 <- function(x1, data, name=data[,1], 
      nicknames=matrix(character(0), 0, 2), ...){
##
de## 1.  Check x1, data, name, nicknames  
##
#  1.1.  x1   
  nx <- length(x1)
  if(nx<1){
    return(list())
  }
  if(!is.character(x1)){
    stop('class(x1) = ', class(x1), '; must be character')
  }
#  1.2.  data
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
#  1.3.  name
  if(!is.character(name)){
    stop('name must be character;  is ', class(name))
  }
#  if(is.numeric(name)){
#    if(length(name)!=1){
#      stop('name is numeric, so its length must be 1;  is ', 
#           length(name))
#    }
#    if(any(name<1)){
#      stop('name is numeric so must be at least 1;  is ', 
#           paste(name[name<1], collapse=', ')) 
#    }
#    if(any(trunc(name)>dimd[2])){
#      stop('name is numeric so must be at most dim(data)[2] = ', 
#        dimd[2], ';  is ', paste(name[name>dimd[2]]) ) 
#    }
#    name <- subNonStandardNames(data[, name], ...)
#  } else if(is.logical(name)){
#    if(length(name) != dimd[2]){
#      stop('name is logical, so its length must match dim(data)[2];', 
#           '  length(name) = ', length(name))
#    }
#    name <- subNonStandardNames(data[, name], ...)
#  } else if(is.character(name)){
  nn <- NROW(name)
  if(nn != dimd[1]){
    stop('NROW(name) = ', nn, ' != nrow(data) = ', 
         nrow(data))
  }
  if(length(dim(name))<2){
    dim(name) <- c(nn, 1)
  }
#      if(NROW(name) != 1){
#        stop("name is character, so its length must be 1;  is ", 
#           length(name))    
#      }
#      sel <- (colnames(data) %in% name) 
#      if(sum(sel) != 1){
#        stop('name is character, so it must be found once in ', 
#           'colnames(data);  found ', sum(sel), ' times')
#      }
#      name <- subNonStandardNames(data[, name], ...)
#    }
#  } else stop('class(name) must be either numeric, logical, or ', 
#              'character;  is ', class(name))
#  1.4  nicknames 
  if(!is.character(nicknames)){
    stop('"nicknames" must be character;  is ', class(nicknames))
  }
  dimn <- dim(nicknames)
  if(length(dimn) != 2){
    ern <- paste0('nicknames must be a matrix; length(dim(nicknames)) = ', 
                  length(dimn))
    stop(ern)
  }
  nrownick <- nrow(nicknames)
  ncol.nick <- ncol(nicknames)
  nrownick2 <- nrownick*ncol.nick
##
## 2.  xsplit <- strsplit(x1, ' ') for first name, middle name, ... 
##  
  xsplit <- strsplit(x1, ' ')
  xlist <- vector(nx, mode='list')
  names(xlist) <- x1
##
## 3.  Process elements of x1 one by one
##
  id <- 1:dimd[1]
  rowi <- function(i, n){
    1+((i-1)%%n)
  }
  for(j in 1:nx){
    xj <- xsplit[[j]]
    jd <- id 
#    dataj <- data 
#    namej <- name
    foundj <- FALSE 
    for(h in seq(along=xj)){
      xj1 <- xj[h]
      matchj2 <- pmatch2(xj1, name[jd, ])[[1]]
      matchj <- rowi(matchj2, nn)
      if(length(matchj)<1){
        xj1s <- subNonStandardNames(xj1)      
        matchj2 <- pmatch2(xj1s, name[jd])[[1]]
        matchj <- rowi(matchj2, nn)
        if(length(matchj)<1) {
#          matchn <- which(nicknames %in% xj1) 
          matchn <- pmatch2(xj1, nicknames)[[1]]
          if(length(matchn)>0){
#            matchn2 <- (1+((matchn+nrownick-1) %% nrownick2))
            matchn2 <- rowi(matchn+nrownick, nrownick2)
            xj1n <- nicknames[matchn2]          
#          matchj. <- (name[jd] %in% xj1n)      
            matchj.2 <- pmatch2(xj1n, name[jd, ])[[1]]
            matchj. <- rowi(matchj.2, nn)
            if(length(matchj.)<1){
              next 
            } else {
              jd <- jd[matchj.]
#            dataj <- data[matchj.,]        
#            namej <- namej[matchj.]
              foundj <- TRUE 
            }
          }
        }
      } else {
        jd <- jd[matchj]
#        dataj <- data[matchj,]
#        namej <- namej[matchj]
        foundj <- TRUE 
      } 
    }
  if(foundj)xlist[[j]] <- jd 
  }  
##
## 4.  Done 
##
  xlist   
}