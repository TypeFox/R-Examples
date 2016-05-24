##################################
#	Internal function called by plot.dtrates(...)
#	Arguments:
#		x,y = coordinates of center of curvature of arc, e.g. (0,0)
#		theta1 = initial theta of arc (radians)
#		theta2 = ending theta of arc (radians)
#		rad = radius of arc
arc <- function(x,y,theta1,theta2,rad,border,...)
{
	noTips <- which((theta2 - theta1) != 0);
	if ((length(theta1)+1)/2 > 1000) {
		steps <- (theta2-theta1)/30; 
		steps <- steps[noTips];
		theta1 <- theta1[noTips];
		theta2 <- theta2[noTips];
		rad <- rad[noTips];
		border <- border[noTips];
		for (i in 1:length(steps))
		{
			xv <- x+rad[i]*cos(seq(theta1[i],theta2[i],steps[i]));
			yv <- y+rad[i]*sin(seq(theta1[i],theta2[i],steps[i]));
			lines(xv,yv,lend=2,col=border[i],...);		
		}
	}
	else { 
	#storing all the coords up front for fast arc plotting, so can be memory intensive.
	#tested on tree with 6670 tips with no problem, but for now only use
	#for trees under 1000 tips
		m <- matrix(NA, nrow=4, ncol=length(noTips));
		m[1,] <- theta2[noTips];
		m[2,] <- theta1[noTips];
		m[3,] <- rad[noTips];
		m[4,] <- border[noTips];
		
		arcsegs <- apply(m, 2, function(z) {
			zz <- as.numeric(z[1:3]);
			inc <- (zz[2] - zz[1])/30
			xv <- zz[3]*cos(seq(zz[1],zz[2],inc));
			yv <- zz[3]*sin(seq(zz[1],zz[2],inc));
			xv <- rep(xv, each=2);
			xv <- xv[-c(1,length(xv))];
			xv <- matrix(xv, ncol=2, byrow=TRUE);
			yv <- rep(yv, each=2);
			yv <- yv[-c(1,length(yv))];
			yv <- matrix(yv, ncol=2, byrow=TRUE);
			data.frame(xv,yv,rep(z[4],nrow(xv)),stringsAsFactors=FALSE);
		});
		arcsegs <- do.call(rbind, arcsegs);
		segments(x+arcsegs[,1], y+arcsegs[,3], x+arcsegs[,2], y+arcsegs[,4], col=arcsegs[,5], lend=2, ...);
	}
}