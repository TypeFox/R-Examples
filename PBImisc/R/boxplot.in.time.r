boxplotInTime <- function (x, xname, additional=T, color = ifelse(additional, "white","lightgrey"), main="", ylim=range(unlist(x),na.rm=T), ..., points = dim(x)[2], at = 1:points) {
    l.obiektow = dim(x)[1]
  col1 = rgb(161, 161, 161, 70, maxColorValue =255)
  col2 = rgb(111, 111, 111, 47, maxColorValue =255)
  col3 = rgb(1, 1, 1, 57, maxColorValue =255)
	plot(range(at),ylim,type="n",xlab="",ylab="",xaxt="n", main=main, ...)
	Axis(at = at, side=1, labels=xname)
	xx = NULL
	for(i in 1:points) {
		if (additional) {
			xx[[i]] = jitter(rep(at[i],times=l.obiektow), amount=0.1)
			}
	}
	for(i in 1:(points-1)) {
		if (additional) {
			segments(xx[[i]], x[,i], xx[[i+1]], x[,i+1], col=col1)
		}                                 
	}
	for(i in 1:points) {
		boxplot(x[,i], add=T, at=at[i],col=color, outline=F)
		if (additional) {
			points(xx[[i]],x[,i],pch=16, col=col3)
			points(xx[[i]],x[,i],pch=16, col=col2,cex=0.5)
			}
	}
}

#boxplot.in.time(kidney[,9:16], colnames(kidney)[9:16], ylim=c(0,110))
