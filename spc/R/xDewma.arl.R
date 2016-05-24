# Computation of EWMA ARLs (drift monitoring)
xDewma.arl <- function(l, c, delta, zr=0, hs=0, sided="one", limits="fix", mode="Gan", m=NULL, q=1, r=40, with0=FALSE) {
  if (l<=0 || l>1) 
    stop("l has to be between 0 and 1")
  if (c<=0) 
    stop("c has to be positive")
  if (zr>c & sided=="one") 
    stop("wrong reflexion border")
  if ( (sided=="two" & abs(hs)>c) | (sided=="one" & (hs<zr | hs>c)) ) 
    stop("wrong headstart")
  if (r<4)
    stop("r is too small")
  ctyp <- pmatch(sided, c("one", "two")) - 1
  if (is.na(ctyp)) 
    stop("invalid ewma type")
  ltyp <- -1 + pmatch(limits, 
          c("fix","vacl","fir","both","Steiner","Knoth","fink","fixW","fixC"))
  if (is.na(ltyp))
    stop("invalid limits type")
  cmode <- pmatch(mode, c("Gan", "Knoth", "Waldmann")) - 1
  if (is.na(cmode))
    stop("invalid algorithm mode")
  if ( is.null(m) ) {
    m <- 0
  } else { if ( m<1 ) stop("m is too small") }
  q <- round(q)
  if (q<1)
    stop("wrong change point position (q)")
  arl <- .C("xDewma_arl",as.integer(ctyp),as.double(l),
            as.double(c),as.double(zr),as.double(hs),
            as.double(delta),as.integer(ltyp),as.integer(m),as.integer(r),
            as.integer(with0),as.integer(cmode),as.integer(q),
            ans=double(length=1),PACKAGE="spc")$ans 
  names(arl) <- "arl"
  return (arl)
}
