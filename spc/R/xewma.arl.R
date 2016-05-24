# Computation of EWMA ARLs (mean monitoring)
xewma.arl <- function(l, c, mu, zr=0, hs=0, sided="one", limits="fix", q=1, r=40) {
  if ( l<=0 | l>2 )
    stop("l has to be between 0 and 2")
  if ( c<=0 )
    warning("usually, c has to be positive")
  if ( zr>c & sided=="one" )
    stop("wrong reflexion border")
  if ( (sided=="two" & abs(hs)>c) | (sided=="one" & (hs<zr | hs>c)) ) 
    warning("unusual headstart")
  if ( r<4 )
    stop("r is too small")
  ctyp <- pmatch(sided, c("one", "two")) - 1
  if ( is.na(ctyp) )
    stop("invalid ewma type")
  ltyp <- -1 + pmatch(limits,
          c("fix", "vacl", "fir", "both", "Steiner", "stat", "fink", "limit", "fixW", "fixC"))
  if ( is.na(ltyp) )
    stop("invalid limits type")
  if ( (sided=="one") & !(limits %in% c("fix", "vacl", "stat", "limit", "fixW")) )
    stop("not supported for one-sided EWMA (not reasonable or not implemented yet")
  q <- round(q)
  if ( q<1 )
    stop("wrong change point position (q)")
  if ( limits=="fix" & q>1 ) {
    arl <- .C("xewma_arl",as.integer(ctyp),as.double(l),
              as.double(c),as.double(zr),as.double(hs),
              as.double(mu),as.integer(ltyp),as.integer(r),as.integer(q),
              ans=double(length=q), PACKAGE="spc")$ans 
  } else {
    arl <- .C("xewma_arl",as.integer(ctyp),as.double(l),
              as.double(c),as.double(zr),as.double(hs),
              as.double(mu),as.integer(ltyp),as.integer(r),as.integer(q),
              ans=double(length=1), PACKAGE="spc")$ans
  }
  names(arl) <- NULL
  return (arl)
}
