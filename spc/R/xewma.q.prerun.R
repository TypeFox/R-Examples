# Computation of EWMA quantiles (mean monitoring) under specified pr-run scenarios
xewma.q.prerun <- function(l, c, mu, p, zr=0, hs=0, sided="two", limits="fix", q=1, size=100, df=NULL, estimated="mu", qm.mu=30, qm.sigma=30, truncate=1e-10, bound=1e-10) {
  if ( l <= 0 | l > 1 )          stop("l (lambda) has to be between 0 and 1")
  if ( c<=0 )                    warning("usually, c has to be positive")
  if ( p <= 0 | p >= 1)          stop("quantile level p must be in (0,1)")
  if ( zr > c & sided == "one")  stop("wrong reflexion border")
  if ( (sided == "two" & abs(hs) > c) | (sided == "one" & ( hs < zr | hs > c )) )
                                 stop("wrong headstart")
  ctyp <- pmatch(sided, c("one", "two")) - 1
  if ( is.na(ctyp) )             stop("invalid ewma type")
  ltyp <- -1 + pmatch(limits, c("fix", "vacl", "fir", "both", "Steiner", "stat"))
  if (is.na(ltyp))               stop("invalid limits type")
  if ( (sided=="one") & !( limits %in% c("fix", "vacl", "stat") ) )
                                 stop("not supported for one-sided EWMA (not reasonable or not implemented yet")
  q <- round(q)
  if ( q<1 )                     stop("wrong change point position (q)")
  if ( size<2 )                  stop("pre run size too small")
  if ( is.null(df) ) df = size - 1
  if ( df<1 )                    stop("degrees of freedom (df) too small")
  emode <- -1 + pmatch(estimated, c("mu", "sigma", "both"))
  if (is.na(emode))              stop("invalid to be estimated type")
  if ( qm.mu<4 )                 stop("qm.mu is too small")
  if ( qm.sigma<4 )              stop("qm.sigma is too small")
  if ( truncate < 0 | truncate >= 0.5 ) stop("wrong value for truncate (should follow 0 < truncate < 0.5)")
  if ( bound < 0 | bound >= 0.001 ) stop("wrong value for bound (should follow 0 < truncate < 0.001)")
  quant <- .C("xewma_q_prerun",
              as.integer(ctyp), as.double(l), as.double(c),
              as.double(p), as.double(zr), as.double(hs),
              as.double(mu), as.integer(ltyp),
              as.integer(q), as.integer(size), as.integer(df), as.integer(emode),
              as.integer(qm.mu), as.integer(qm.sigma), as.double(truncate), as.double(bound),
              ans=double(length=1),PACKAGE="spc")$ans
  names(quant) <- "q"
  quant
}
