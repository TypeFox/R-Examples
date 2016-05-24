checkNames <- function(x, warn=0, unique=TRUE, 
       avoid=character(0), 
       message0=head(deparse(substitute(x), 25), 2), 
       ...){
##
## 1.  names(x)
##
  namex <- names(x)
  na <- length(avoid)
  if((na%%2)!=0){
    stop('avoid = ', deparse(substitute(avoid), 25), 
         '\nmust have length = a nonnegative even integer;  is ',
         na)    
  }
  ja <- seq_len(na/2)
#  
  makeNames <- function(y, u){
    n.y <- make.names(y, u)
    for(j in ja){
      m <- which((regexpr(avoid[2*j-1], n.y) > 0) 
                  & (y != n.y) )
      if(length(m)>0){
        n.y[m] <- sub(avoid[2*j-1], avoid[2*j], 
                      n.y[m])
        m2 <- grep(avoid[2*j-1], n.y)
        if(length(m2)>0){
          msg <- paste0('Cannot fix duplicate name ', 
            'for x[', m2, '];  names(x)[', m2, ']') 
          if(nchar(y[m2])<1){
            msg2 <- paste0(msg, ' was blank, and ', 
                'make.names could not fix it with ',
                'avoid = ', paste(avoid, collapse = '; '))            
          } else {
            msg2 <- paste0(msg, ' = ', y[m2], 
                ' converted with avoide to ', n.y[m2])
          }
          stop(msg2)
        }
        next      
      }
    }
    n.y 
  }
##
## 2.  NULL?
##
  if(is.null(namex)){
    msg0 <- paste0(paste0(message0, collapse=":"), 
        ":  names = NULL; returning make.names(character(", 
        "length(x))), ", 
        unique, ")")
    if(warn>(-1)){
      op <- options(warn=warn)
      on.exit(options(op))
      warning(msg0)
    }
    nmx0 <- makeNames(character(length(x)), unique)
#
    attr(nmx0, 'message') <- msg0
    return(nmx0) 
  }
##
## 3.  unique?
##  
  if(unique){
    tabNames <- table(namex)
    if(length(tabNames) < length(x)){
      tN2 <- tabNames[tabNames>1][1]
      msg <- paste0(paste(message0, collapse=":") , 
          ":  names(x) not unique; ", 
          names(tN2), " occurs ", tN2, 
          ' times; returning make.names(names(x), TRUE)')
      if(warn>(-1)){
        op <- options(warn=warn)
        on.exit(options(op))
        warning(msg)
      }
      nameu <- makeNames(namex, TRUE)
      attr(nmeu, 'message') <- msg
      return(nmeu)   
    }
  }
  namex 
}