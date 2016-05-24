lval <-
function (x, na.rm = TRUE) 
{
  xna <- is.na(x)
  if (na.rm) 
    x <- x[!xna]
  else if (any(xna)) 
    return(rep(NA, 5))
  x <- sort(x)
  n <- length(x)
  cpos <- n
  depth <- c()
  while (cpos > 1) {
    cpos <- (floor(cpos) + 1)/2
    if (cpos != 1.5) 
      depth <- c(depth, cpos)
  }
  lo <- (x[floor(depth)] + x[ceiling(depth)])/2
  hi <- (x[floor(n + 1 - depth)] + x[ceiling(n + 1 - depth)])/2
  mids <- (lo + hi)/2
  spreads <- hi - lo
  out = data.frame(depth, lo, hi, mids, spreads)
  labels = c("M", "H", "E", "D", "C", "B", "A", "Z", "Y", "X", 
             rep("", 1000))
  row.names(out) = labels[1:length(dimnames(out)[[1]])]
  return(out)
}
