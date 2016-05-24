mkDummy <- function(x) {

  xN <- deparse(substitute(x))
  if (!inherits(x,"factor")) stop("x should be a factor")
  lev <- levels(x)

  result <- lapply(lev,function(l) as.integer(x==l))
  names(result) <- paste(xN,".",lev,sep="")
  result

}
