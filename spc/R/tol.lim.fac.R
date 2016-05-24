# Computation of 2-sided tolerance limits factors
tol.lim.fac <- function(n,p,a,mode="WW",m=30) { 
  if (n<2) 
    stop("n has to be larger than 1")
  if (p<=0 | p>=1) 
    stop("p has to be in (0,1)")
  if (a<=0 | a>=1)
    stop("a has to be in (0,1)")
  mtype <- pmatch(mode, c("WW", "exact")) - 1
  if (is.na(mtype))
    stop("invalid mode type")
  if (m<10)
    stop("m has to be at least 10")
  tlf <- .C("tol_lim_fac",as.integer(n),as.double(p),
             as.double(a),as.integer(mtype),as.integer(m),
             ans=double(length=1),PACKAGE="spc")$ans 
  names(tlf) <- "k"
  return (tlf)
}
