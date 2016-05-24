"aq.plot" <-
function(x, delta=qchisq(0.975, df=ncol(x)), quan=1/2, alpha=0.05) {

  if(is.vector(x) == TRUE || ncol(x) == 1) { stop("x must be at least two-dimensional") }

  covr <- covMcd(x, alpha=quan)
  dist <- mahalanobis(x, center=covr$center, cov=covr$cov)
  s <- sort(dist, index=TRUE)

  z <- x
  if(ncol(x) > 2) {
	p <- princomp(x,covmat=covr)
	z <- p$scores[,1:2]
	sdprop <- (p$sd[1]+p$sd[2])/sum(p$sd)
	cat("Projection to the first and second robust principal components.\n")
	cat("Proportion of total variation (explained variance): ")
	cat(sdprop)
	cat("\n")
  }
	
    par(mfrow=c(2,2), mai=c(0.8,0.6,0.2,0.2), mgp=c(2.4,1,0))
    plot(z, col=3, type="n", xlab="", ylab="")
    text(z, dimnames(as.data.frame(z))[[1]], col=3, cex=0.8)

  plot(s$x, (1:length(dist))/length(dist), col=3, xlab="Ordered squared robust distance", ylab="Cumulative probability", type="n")
  text(s$x, (1:length(dist))/length(dist), as.character(s$ix), col=3, cex=0.8)
  t <- seq(0,max(dist), by=0.01)
  lines(t, pchisq(t, df=ncol(x)), col=6)

  abline(v=delta, col=5)
  text(x=delta, y=0.4, paste(100*(pchisq(delta,df=ncol(x))),"% Quantile",sep=""), col=5, pos=2, srt=90, cex=0.8)

  xarw <- arw(x, covr$center, covr$cov, alpha=alpha)
  if (xarw$cn<Inf){
    abline(v=xarw$cn, col=4)
    text(x=xarw$cn, y=0.4, "Adjusted Quantile", col=4, pos=4, srt=90, cex=0.8)
  }
    plot(z, col=3, type="n", main=paste("Outliers based on ",100*(pchisq(delta,df=ncol(x))),"% quantile",sep=""), xlab="", ylab="")
    if (any(dist>delta)){
	text(z[dist>delta, 1], z[dist>delta, 2],dimnames(as.data.frame(x)[dist>delta,])[[1]],col = 2, cex = 0.8)
    }
    if (any(dist<=delta)){
	text(z[dist<=delta, 1], z[dist<=delta, 2], dimnames(as.data.frame(x)[dist<=delta,])[[1]],col = 3, cex = 0.8)
    }
    plot(z, col=3, type="n", main="Outliers based on adjusted quantile", xlab="", ylab="")
  if (xarw$cn<Inf){
    text(z[dist>xarw$cn, 1], z[dist>xarw$cn, 2], dimnames(as.data.frame(x)[dist>xarw$cn,])[[1]],col = 2, cex = 0.8)
  }
    text(z[dist<=xarw$cn, 1], z[dist<=xarw$cn, 2], dimnames(as.data.frame(x)[dist<=xarw$cn,])[[1]],col = 3, cex = 0.8)
   o <- ( sqrt(dist) > max(sqrt(xarw$cn), sqrt(qchisq(0.975, dim(x)[2]) ) ) )
   list(outliers = o)
}

