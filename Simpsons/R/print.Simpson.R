print.Simpson <-
function(x,...)
{
  Res <- x$alldata
  names(Res)[1] <- x$namex
  names(Res)[2] <- x$namey
  return(Res)
}
