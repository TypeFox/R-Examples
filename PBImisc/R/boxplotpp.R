`boxplotpp` <-
function(x, xname=seq(1:ncol(x)), utitle="", addLines=TRUE, color = ifelse(addLines, "white","lightgrey"), ...) {
  l.punktow = ncol(x)
  l.obiektow = nrow(x)
	plot(c(0.5,l.punktow+0.5),range(x),type="n",xlab="",ylab="",xaxt="n", main=utitle, ...)
	Axis(at = 1:l.punktow, side=1, labels=xname)
	xx = NULL
	for(i in 1:l.punktow) {
		if (addLines) {
			xx[[i]] = jitter(rep(i,times=l.obiektow), amount=0.0625)
			}
	}
	for(i in 1:(l.punktow-1)) {
		if (addLines) {
			segments(xx[[i]], x[,i], xx[[i+1]], x[,i+1], col="lightgrey")
		}
	}
	for(i in 1:l.punktow) {
		boxplot(x[,i], add=T, at=i,col=color)
		if (addLines) {
			points(xx[[i]],x[,i],pch=16, col="black")
			points(xx[[i]],x[,i],pch=16, col="lightgrey",cex=0.5)
			}
	}
}

