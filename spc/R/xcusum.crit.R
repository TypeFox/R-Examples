# Computation of CUSUM decision limits for given ARL (mean monitoring)
xcusum.crit <- function(k, L0, mu0 = 0, hs = 0, sided = "one", r = 30) {
  if (k<0)
    stop("k has to be non-negative")
  if (L0<1) 
    stop("L0 is too small")
  if (hs<0) 
    stop("wrong headstart")
  if (r<4) 
    stop("r is too small")
  ctyp <- pmatch(sided, c("one", "two", "Crosier")) - 1
  if (is.na(ctyp)) 
    stop("invalid cusum type")
  h <- .C("xcusum_crit",as.integer(ctyp),as.double(k),
          as.double(L0),as.double(hs),as.double(mu0),as.integer(r),
          ans=double(length=1),PACKAGE="spc")$ans 
  names(h) <- "h"
  return (h)
}

