hr05CutoffMvnormal <- 
# 
# corrected cutoff values, degrees of freedom, etc.
# from Hardin Rocke 2005 paper
# assumes MCD subset comes from MV normal
#
# Author: Christopher G. Green
# Date: 2010-02-05
#
function( n.obs, p.dim, 
  mcd.alpha=max.bdp.mcd.alpha(n.obs, p.dim), signif.alpha=0.05, 
  method=c("GM14","HR05"), use.consistency.correction=FALSE ) 
{

  method          <- match.arg(method)

  ch99             <- ch99AsymptoticDF( n.obs, p.dim, mcd.alpha )

  # use hardin and rocke 2005 results or my simulated versions
  # to get predicted degrees of freedom value

  m.pred           <- hr05AdjustedDF(n.obs=n.obs, p.dim=p.dim, 
    mcd.alpha=mcd.alpha, m.asy=ch99$m.hat.asy, method=method)

  cutoff.pred      <- hr05CriticalValue(m.pred        ,p.dim,signif.alpha)
  cutoff.asy       <- hr05CriticalValue(ch99$m.hat.asy,p.dim,signif.alpha)

  if ( use.consistency.correction ) {
    cutoff.pred <- cutoff.pred * ch99$c.alpha
    cutoff.asy  <- cutoff.asy  * ch99$c.alpha
  }

  list( cutoff.pred = cutoff.pred, cutoff.asy = cutoff.asy, 
    c.alpha = ch99$c.alpha, m.asy = ch99$m.hat.asy, m.pred = m.pred, 
    n = n.obs, p = p.dim )
}
