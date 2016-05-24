# Computation of GRSR (Girshick, Rubin, Shiryaev, Roberts) steady-state ARLs (mean monitoring)
xgrsr.ad <- function(k, g, mu1, mu0=0, zr=0, sided="one", MPT=FALSE, r=30) {
  if (k<0) 
    stop("k has to be non-negative")
  if (g<0)
    stop("g has to be positive")
  if (r<4) 
    stop("r is too small")
  ctyp <- pmatch(sided, c("one", "two")) - 1
  if (is.na(ctyp)) 
    stop("invalid grsr type")
  ad <- .C("xgrsr_ad",as.integer(ctyp),as.double(k),
           as.double(g),as.double(mu0),as.double(mu1),as.double(zr),as.integer(r),as.integer(MPT),
           ans=double(length=1),PACKAGE="spc")$ans 
  names(ad) <- "ad"
  return(ad)
}
