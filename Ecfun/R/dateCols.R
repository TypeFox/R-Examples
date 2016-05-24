dateCols <- function(col.names, YMD=c('Year', 'Month', 'Day')){
##
## 1.  colNms 
##
  colNms <- colnames(col.names)
  if(is.null(colNms))
    colNms <- col.names   
##
## 2.  ymd <- grep(YMD, colNms)
##
  if(length(YMD)!=3)
    stop('YMD must have length 3;  is ', length(YMD))
#  get subscripts of desired elements of colNms s
  ymd <- lapply(YMD, grep, x=colNms)
##
## 3.  gpNames <- sub(YMD, '', ymd) 
##  
  gpNames <- vector(mode='list', 3)
  names(gpNames) <- YMD
  for(i in 1:3){
    gpNames[[i]] <- sub(YMD[i], '', colNms[ymd[[i]]])
  }
##
## 4.  Try to match gpNames[[1]] with [[2]] and [[3]]
##
  gps <- gpNames[[1]]
  ng <- length(gps)
#  4.1.  check lengths 
  n2 <- length(gpNames[[2]])
  if(n2 != ng)
    warning('number of matches for Year = ', ng, 
            ' != number of matches for Month = ', n2) 
  n3 <- length(gpNames[[3]])
  if(n3 != ng)
    warning('number of matches for Year = ', ng, 
            ' != number of matches for Day = ', n3) 
#
  out <- vector(mode='list', ng)
  names(out) <- gps
  for(j in 1:ng){
    iy <- which(gps==gps[[j]])
    oj <- ymd[[1]][iy]
#    
    im <- which(gpNames[[2]]==gps[[j]])
    if(length(im)<1){
      warning('No Month found to match Year ', gps[[j]])
      next 
    } else if(length(im)>1){
      warning('multiple Months found to match Year ', gps[[j]], 
              ';  ignoring')
      next 
    }
    oj2 <- c(oj, ymd[[2]][im])  
#
    id <- which(gpNames[[3]]==gps[[j]])
    if(length(id)<1){
      warning('no Day found to match Year ', gps[[j]])
      next
    } else if(length(id)>1){
      warning('multiple Days found to match Year ', gps[[j]], 
              ';  ignorning')
      next
    }
    oj3 <- c(oj2, ymd[[3]][id])
    names(oj3) <- YMD
    out[[j]] <- oj3
  }
##
## 5.  Return a list of integer vectors of length 3 for each 
##  triple found.  
##  
  out 
}
