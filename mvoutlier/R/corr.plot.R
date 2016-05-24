"corr.plot" <-
function(x, y, quan=1/2, alpha=0.025, ...) {

  #library(rrcov)  
  x <- as.matrix(cbind(x,y))

  covr <- covMcd(x, cor=TRUE, alpha=quan)
  cov.svd <- svd(cov(x), nv = 0)
  covr.svd <- svd(covr$cov, nv = 0)
  r <- cov.svd[["u"]] %*% diag(sqrt(cov.svd[["d"]]))
  rr <- covr.svd[["u"]] %*% diag(sqrt(covr.svd[["d"]]))
  
  e <- cbind(cos(c(0:100)/100 * 2 * pi) * sqrt(qchisq(1-alpha,2)), sin(c(0:100)/100 * 2 * pi) * sqrt(qchisq(1-alpha,2)))
  tt <- t(r %*% t(e)) + rep(1, 101) %o% apply(x, 2, mean)
  ttr <- t(rr %*% t(e)) + rep(1, 101) %o% covr$center
  plot(x, xlim=c(min(c(x[, 1], tt[, 1],ttr[,1])), max(c(x[, 1], tt[, 1],ttr[,1]))), ylim=c(min(c(x[, 2], tt[, 2],ttr[,2])), max(c(x[, 2], tt[, 2],ttr[,2]))), ...)
  title(main=list(paste("Classical cor =",round(cor(x)[1,2],2),"                                      "),col=4))
  title(main=list(paste("                                      Robust cor =",round(covr$cor[1,2],2)),col=2))

  lines(tt[, 1], tt[, 2], type = "l",col=4,lty=3)
  lines(ttr[, 1], ttr[, 2], type = "l",col=2)

  ret <- list(cor.cla = cor(x)[1,2], cor.rob=covr$cor[1,2])
  ret
}

