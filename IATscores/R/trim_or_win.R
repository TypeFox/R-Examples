trim_or_win <- function (x, type = "trm", tr = .1, na.rm = FALSE)
{
  # trim_or_win is a funciton that trims or winsorizes data according
  # to parameters typemean and trmean
  if(!type %in% c("std", "trm", "wns", "inv"))
    return("ERROR: invalid parameter type")
  
  netx <- x[!is.na(x)]
  
  if (type == "std")
  {
    newx <- netx
  }
  if (type == "trm")
  {
    newx <- trimval(netx, tr = tr)
  } else if (type == "wns")
  {
    newx<-winval (netx, tr = tr)
  } else if (type == "inv")
  {
    newx <- inv_trimval(netx, tr = tr)
  }
  x[!is.na(x)] <- newx
  if(na.rm == TRUE) x <- x[!is.na(x)]
  x
}