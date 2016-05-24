# Computation of EWMA quantiles (mean monitoring, t distributed data)
xtewma.q <- function(l, c, df, mu, alpha, zr=0, hs=0, sided="two", limits="fix", mode="tan", q=1, r=40) {
  if ( l <= 0 | l > 1 )          warning("l is typically between 0 and 1 -- you should really know what you do")
  if ( c<=0 )                    warning("usually, c has to be positive")
  if ( df < 1 )                  stop("df must be greater or equal to 1")
  if ( alpha <= 0 | alpha >= 1)  stop("quantile level alpha must be in (0,1)")
  if ( zr > c & sided == "one")  stop("wrong reflexion border")
  if ( (sided == "two" & abs(hs) > c) | (sided == "one" & ( hs < zr | hs > c )) )
                                 warning("unusual headstart")
  ctyp <- pmatch(sided, c("one", "two")) - 1
  if ( is.na(ctyp) )             stop("invalid ewma type")
  ltyp <- -1 + pmatch(limits, c("fix", "vacl"))
  if (is.na(ltyp))               stop("invalid limits type")
  ntyp <- -1 + pmatch(mode, c("identity", "sin", "sinh", "tan"))
  if ( is.na(ntyp) )             stop("substitution type not provided (yet)") 
  if ( r < 4 )                   stop("r is too small")
  q <- round(q)
  if ( q<1 )                     stop("wrong change point position (q)")
  quant <- .C("xtewma_q",
              as.integer(ctyp), as.double(l), as.double(c), as.double(alpha), as.double(zr),
              as.double(hs), as.integer(df), as.double(mu),
              as.integer(ltyp), as.integer(r), as.integer(ntyp), as.integer(q),
              ans=double(length=1),PACKAGE="spc")$ans
  names(quant) <- "q"
  quant
}
