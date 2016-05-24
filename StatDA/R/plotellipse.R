"plotellipse" <-
function(x.loc, x.cov, perc=0.98, col=NULL, lty=NULL)
{
# plots perc % tolerance ellipses with certain location and covariance
# perc should be a (vector of) number(s) between 0 and 1

  if (ncol(x.cov)!=2) stop("Dimension of data must be 2!")

  if (is.null(col)) col <- rep(1,length(perc))
  if (is.null(lty)) lty <- rep(1,length(perc))

  cov.svd <- svd(x.cov, nv = 0)
  r <- cov.svd[["u"]] %*% diag(sqrt(cov.svd[["d"]]))

  m <- 100

  for (i in 1:length(perc)){
    alphamd=sqrt(qchisq(perc[i], 2))
    e1md <- cos(c(0:m)/m * 2 * pi) * alphamd
    e2md <- sin(c(0:m)/m * 2 * pi) * alphamd
    emd <- cbind(e1md, e2md)
    ttmd <- t(r %*% t(emd)) + rep(1, m + 1) %o% x.loc
    lines(ttmd[, 1],ttmd[, 2], col=col[i], lty=lty[i])
  }
}

