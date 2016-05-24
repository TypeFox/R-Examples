strsplit1 <- function(x, split=',', Quote='"', ...){
##
## 1.  spl1 <- regexpr(split, x, ...)
##
  spl1 <- regexpr(split, x, ...)
##  
## 2.  Qt1 <- regexpr(Quote, x, ...)
##
  Qt1 <- regexpr(Quote, x, ...)
##
## 3.  Qte <- (Qt1<spl1)
##
  Qte <- ((0<Qt1) & (Qt1<spl1))
  Spl1 <- spl1 
  if(any(Qte)){
#   quote followed by a split     
    Spl1[Qte] <- (-1)
    Qt1. <- Qt1[Qte]
    xQt1 <- substr(x[Qte], Qt1.+1, nchar(x[Qte]))
    Qt2 <- regexpr(Quote, xQt1, ...)
    if(any(Qt2>0)){
#     matching quote       
      xQt2 <- xQt1[Qt2>0]
      Qt2. <- Qt2[Qt2>0]
      xQt2. <- substring(xQt2, Qt2.+1, nchar(xQt2))
      spl.1 <- regexpr(split, xQt2.)
      if(any(spl.1>0)){
#        Need this in case 2 commas, 
#        one before and one after the close quote:           
        Spl.1 <- (Qt1.[Qt2>0]+Qt2.+spl.1)[spl.1>0]
        Spl1[Qte][Qt2>0][spl.1>0] <- Spl.1 
      } 
    }
  }
##
## 4.  x1  
##
  x1 <- x
  x1[Spl1==1] <- ''
  x1[Spl1>1] <- substr(x[Spl1>1], 1, Spl1[Spl1>1]-1)
##
## 5.  x2 
##
  x2 <- x
  x2[Spl1<0] <- ''
  x2[Spl1>0] <- substr(x[Spl1>0], Spl1[Spl1>0]+1, nchar(x[Spl1>0]))
##
## 6.  Done 
##
  list(x1, x2)  
}
