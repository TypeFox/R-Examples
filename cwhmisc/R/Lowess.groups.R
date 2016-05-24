lowess.bygroup <- function(x, y, group, span=2/3, col=seq_along(x), lty=seq_along(x)) {
  ind <- sort.list(x)
  x <- x[ind]; y <- y[ind]
  lc <- length(col); lt <- length(lty)
  ly <- length(y)
  if (lc == 1) col <- rep(col,ly)
  if (lt == 1) col <- rep(lty,ly)
  for (ii in unique(group)) {
    iig <- ii == group
    xi  <- x[iig]; yi <- y[iig]
    lines(lowess(xi,yi,f=span),col=col[ii],lty=lty[ii])
  }
}
