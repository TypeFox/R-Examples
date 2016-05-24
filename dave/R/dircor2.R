dircor2<- function(veg,x.axis,y.axis,step=5) {
	D.veg <- as.dist((1 - cor(t(veg)))/2)    # Correlation, transformed to distance
	coor <- cbind(x.axis,y.axis)                  
	steps <- seq(0,180,step)
	ll <- length(steps)
	rout <- rep(0,ll)                       # mean r
	rlow <- rep(0,ll)                       # lower limit of r
	rupp <- rep(0,ll)                       # upper limit of r
	E <- dist(coor,method="euclidean")
	M <- dist(coor[,1],method="manhattan")
	for(i in 1:ll){   
		phid <- 2*pi*(steps[i]/360)
		Alpha <- acos(M/E)
		Beta <- Alpha-phid
		P <- E*cos(Beta)
		mstat <- mantel(D.veg,P)
		rout[i] <- mstat$statistic
		rlow[i] <- rout[i]+sort(mstat$perm)[50]  # heuristic confidence limit
		rupp[i] <- rout[i]+sort(mstat$perm)[950] # heuristic confidence limit
	}
	plot(c(steps,steps),c(rlow,rupp),xlab="Direction, deg",ylab="Mantel correlation",type="n")
	points(steps,rout,pch=1,cex=1.0)
	lines(steps,rout,lwd=0.5)
	lines(steps,rlow,lwd=0.3)
	lines(steps,rupp,lwd=0.3)
	abline(v=c(50,100,150),lwd=1,col="gray")
	abline(h=0,lwd=1,col="gray")
# output to list dircor
	o.dircor<- list(method="directional mantel correlation",steps=steps,mean.correlation=rout,lower.limit=rlow,upper.limit=rupp)
}
