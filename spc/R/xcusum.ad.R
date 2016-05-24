# Computation of CUSUM steady-state ARLs (mean monitoring)
xcusum.ad <- function(k, h, mu1, mu0 = 0, sided = "one", r = 30) {
  if (k<0) 
    stop("k has to be non-negative")
  if (h<=0) 
    stop("h has to be positive")
  if (r<4) 
    stop("r is too small")
  if (r>30 & r<=50 & sided=="two") 
    warning("computation needs some time")
  if (r>50 & sided=="two") 
    warning("ought to be restricted to very fast CPUs")
  ctyp <- pmatch(sided, c("one", "two", "Crosier")) - 1
  if (is.na(ctyp)) 
    stop("invalid cusum type")
  ad <- .C("xcusum_ad",as.integer(ctyp),as.double(k),
           as.double(h),as.double(mu0),as.double(mu1),as.integer(r),
           ans=double(length=1),PACKAGE="spc")$ans 
  names(ad) <- "ad"
  return(ad)
}
