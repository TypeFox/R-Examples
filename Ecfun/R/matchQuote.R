matchQuote <- function(x, Quote='"', sep=' ', 
                       maxChars2append=2, ...){
##
## 1.  nch, gsub, odd number of substitutions?  
##
  nchx <- nchar(x)
  xq <- gsub(Quote, '', x, ...)
  nchq <- nchar(xq)
  dn <- (nchx-nchq)
  osub <- which((dn%%2) > 0) 
  xo <- x
#  attr(xo, 'unmatchedQuotes') <- osub 
##
## 2.  blank lines to drop?   
##
  nx <- length(x)
# os1 = (osub+1), dropping any osub==nx 
  osx <- osub[osub<nx]
  if(length(osx)<1){
#   lines with an unmatched Quote before the last line  
    attr(xo, 'unmatchedQuotes') <- osub     
    attr(xo, 'blankLinesDropped') <- integer(0)
    attr(xo, 'quoteLinesAppended') <- integer(0)
    attr(xo, 'ncharAppended') <- integer(0)
    return(xo)
  }
  os1 <- (osx+1)
  xnb <- gsub(' ', '', x[os1])
  bld <- os1[nchar(xnb)<1] 
# attr(xo, 'blankLinesDropped') <- bld 
##
## 3.  next line quote?
##
#  3.1.  find 
  os2 <- os1
  os2[nchar(xnb)<1] <- (os1[nchar(xnb)<1]+1)
  unmatch2.0 <- which(os2 %in% osx)
  onmatched <- os2[unmatch2.0]
  unmatch2 <- which(nchx[onmatched] <= maxChars2append) 
  Unmatch2 <- onmatched[unmatch2]
#  Os2 <- os2[Unmatch2]
# attr(xo, 'quoteLinesAppended') <- Os2 <- Unmatch2  
# attr(xo, 'ncharsAppended') <- nchx[Unmatch2]
##
## 4.  append 
##
  if(length(Unmatch2)>0){
    xo[osx[unmatch2]] <- paste(x[osx[unmatch2]], 
                               x[Unmatch2], sep=sep)
  }
##
## 5.  Drop
##
#  Drop <- c(bld, Os2)
  Drop <- c(bld, Unmatch2)
  if(length(Drop)>0){ 
    Xo <- xo[-Drop]
  } else Xo <- xo 
##
## 6.  Done 
##
  attr(Xo, 'unmatchedQuotes') <- osub 
  attr(Xo, 'blankLinesDropped') <- bld 
  attr(Xo, 'quoteLinesAppended') <- Unmatch2  
  attr(Xo, 'ncharsAppended') <- nchx[Unmatch2]
  Xo   
}
