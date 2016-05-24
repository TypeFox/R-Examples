as.POSIXct1970 <- function(x, tz="GMT", ...){
#
  if(!is.numeric(x)){
    Px <- try(as.POSIXct(x, tz=tz, ...))
    if(class(Px)[1]=='try-error'){
      nx <- length(x)
      Px <- rep(as.POSIXct1970(0), nx)
    }
    return(Px)
  }
#
  o1970 <- strptime('1970-01-01', '%Y-%m-%d', tz=tz)
  o1970. <- as.POSIXct(o1970, tz=tz)
#
  as.POSIXct(x, tz=tz, origin=o1970., ...)
}
