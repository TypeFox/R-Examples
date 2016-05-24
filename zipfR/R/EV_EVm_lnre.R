## generic stubs for EV and EVm (when not provided by the respective LNRE model
EV.lnre <- function (obj, N=NA, ...)
{
  if (! inherits(obj, "lnre")) stop("argument must be object of class 'lnre'")
  stop("EV() method not implemented for ",obj$name," LNRE model (lnre.",obj$type,")")
}

EVm.lnre <- function (obj, m, N=NA, ...)
{
  if (! inherits(obj, "lnre")) stop("argument must be object of class 'lnre'")
  stop("EVm() method not implemented for ",obj$name," LNRE model (lnre.",obj$type,")")
}
