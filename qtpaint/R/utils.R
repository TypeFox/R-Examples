recycleVector <- function(x, len) {
  xlen <- length(x)
  if (xlen == len)
    return(x)
  if (!len)
    return(rep(x, 0))
  if (xlen > len) {
    warning("the length of 'x' is already longer than 'len'")
    return(x)
  }
  if (!xlen)
    stop("cannot recycle a zero length vector")
  times <- as.integer(len %/% xlen)
  if (times * xlen != len) {
    warning("longer argument not a multiple of shorter argument length")
    return(rep(x, length = len))
  }
  rep.int(x, times) # 100 X faster when length(x) is 1
}
