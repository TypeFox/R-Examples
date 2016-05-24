fnSubset <- function(x, fnFull, xFixed, xFull=c(x, xFixed), ...){
##
## 1.  Confirm length(x)+length(xFixed) = length(xFull)
##
  nx <- length(x)
  nFixed <- length(xFixed)
  nFull <- length(xFull)
  if((nx+nFixed) != nFull)
    stop("length(x)+length(xFixed) != length(xFull):  ",
         nx, " + ", nFixed, " != ", nFull)
##
## 2.  names(xFull)?
##
#  2.1.  is.null(names(xFull))  
  if(is.null(names(xFull)))
    return(fnFull(c(x, xFixed), ...))
#  2.2.  xFull[names(xFixed)] <- xFixed, ... 
  {
    if(is.null(names(xFixed))){
      if(is.null(names(x)))
        xFull <- c(x, xFixed)
      else {
        x. <- (names(xFull) %in% names(x))
        if(sum(x.) != nx){
          print(x)
          print(xFull) 
          stop("x has names not in xFull.")
        }
        xFull[names(x)] <- x
        xFull[!x.] <- xFixed 
      }
    }
    else {
      Fixed <- (names(xFull) %in% names(xFixed))
      if(sum(Fixed) != nFixed){
        print(xFixed)
        print(xFull)
        stop("xFixed has names not in xFull.")
      }
      xFull[names(xFixed)] <- xFixed
      {
        if(is.null(names(x))) xFull[!Fixed] <- x
        else {
          x. <- (names(xFull) %in% names(x))
          if(sum(x.) != nx){
            print(x)
            print(xFull) 
            stop("x has names not in xFull.")
          }
          xFull[names(x)] <- x
        }
      }
    }
  }
##
## 3.  fnFull(...)
##
  fnFull(xFull, ...)
}
