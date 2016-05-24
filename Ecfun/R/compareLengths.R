compareLengths <- function(x, y, 
          name.x=deparse(substitute(x), width.cutoff, nlines=1, ...), 
          name.y=deparse(substitute(y), width.cutoff, nlines=1, ...), 
          message0='', compFun=c('NROW', 'length'), 
          action=c(compatible='', incompatible='warning'), 
          length0=c('compatible', 'incompatible', 'stop'), 
          width.cutoff=20, ...){
##
## 1.  nchar(name.x, name.y)?
## 
  if((nchar(name.x)<1) || (nchar(name.y)<1)){
    message0 <- paste0(message0, 'in compareLengths:')
  }
  if(nchar(name.x)<1) name.x <- 'x'
  if(nchar(name.y)<1) name.y <- 'y'  
##
## 2.  lenx, leny
##
  comp <- match.arg(compFun)
  lenx <- do.call(comp, list(x))
  if(length(lenx)!=1){
    stop(message0, ' compFun[ = ', comp, '](', name.x, ') has length ', 
         lenx, '; must be 1.')    
  }
  if(!is.numeric(lenx)){
    stop(message0, ' compFun[ = ', comp, '](', name.x, 
         ') is not numeric;  class = ', class(lenx)) 
  }
#
  leny <- do.call(comp, list(y))
  if(length(leny)!=1){
    stop(message0, ' compFun[ = ', comp, '](', name.y, ') has length ', 
       leny, '; must be 1.')    
  }
  if(!is.numeric(leny)){
    stop(message0, ' compFun[ = ', comp, '](', name.y, 
         ') is not numeric;  class = ', class(leny)) 
  }
  len <- c(lenx, leny)
##
## 3.  lenx==leny?
##  
  if(lenx==leny)return(c('equal', ''))
##
## 4.  lenx=0 or leny=0?    
##
  act <- match.arg(action)
  o <- order(len)
  nam <- c(name.x, name.y)
  res <- (len[o[2]] %% len[o[1]])
  if(is.na(res)){
    ms0 <- paste0(message0, ' length(', nam[o[1]], 
                 ') = 0')
    if(length0[1] == 'compatible'){
      Ms0 <- c('compatible', ms0)
      if(nchar(action[1])<1){ 
        return(Ms0)
      } else do.call(action[1], list(Ms0))
    } else {
      if(length0[1]=='stop'){
        stop(ms0)
      }
      Ms0 <- c('incompatible', ms0)
      if(nchar(action[2])<1){
        return(Ms0)
      } else do.call(action[2], list(Ms0))
    }
  }
##
## 5.  compatible? 
##
  if(res==0){
    rat <- (len[o[2]] %/% len[o[1]])
    msc <- paste0(message0, ' length(', nam[o[2]], ') = ', 
        len[o[2]], ' is ', rat, ' times length(', 
        nam[o[1]], ') = ', len[o[1]])
    Msc <- c('compatible', msc)
    if(nchar(action[1])<1){ 
      return(Msc)
    } else do.call(action[1], list(Msc))
  }  
##
## 6.  incompatible   
##
  msi <- paste0(message0, ' length(', nam[o[2]], ') = ', 
        len[o[2]], ' is not a multiple of length(', 
        nam[o[1]], ') = ', len[o[1]])
  Msi <- c('incompatible', msi)
  if(nchar(action[2])<1){ 
    return(Msi)
  } else { 
    do.call(action[2], list(Msi[1], ': ', Msi[2]))
  }
  Msi
}
