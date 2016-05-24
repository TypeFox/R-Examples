`plotCuspBifurcation` <-
function(object, xlim=a+c(-0.3,0.3), ylim=b+c(-.1,0.1), xlab = expression(alpha),
	ylab=expression(beta), hue=0.5+0.25*tanh(object$y), 
	col=hsv(h=hue,s=1,alpha=.4), cex.xlab=1.55, cex.ylab=cex.xlab, 
	axes=TRUE, box = TRUE, add=FALSE, bifurcation.set.fill=gray(.8), cex.scale = 15, 
	cex = (cex.scale/log(NROW(ab)))*dens/max(dens), pch=20 ){
	obj <- object
	ab <- obj$linear.predictors
	a <- c(-1,1) * max(abs(range(ab[,'alpha'])))
	b <- c(-1,1) * max(abs(range(ab[,'beta'])))
	if(a[2]>b[2]) {b <- a} else {a <- b}
	if(!add){
		plot.new();
		plot.window(xlim,ylim);
		if(axes) {
			axis(1); 
			axis(2); 
		}
		if(box){
			box();
		}
		mtext(xlab,1,line=2,cex=cex.xlab);
		mtext(ylab,2,line=3,cex=cex.ylab);
		bif <- cusp.bifset(seq_range(b+c(-.7,.7)))
		polygon(c(bif[,2],rev(bif[,3])),c(bif[,1],rev(bif[,1])),col=bifurcation.set.fill)
		abline(v=0, lty=3, col=gray(.3))
		abline(h=0, lty=3, col=gray(.3))
	} 

	density2D = Vectorize((function(x,y,bw=.9)sum(exp(-0.5*colSums((t(ab)-c(x,y))^2/bw^2)))/sqrt(2*pi*bw^2)/NROW(ab)))
	dens = density2D(ab[,'alpha'],ab[,'beta'])
	#hcl(240,l=100-100*dens/max(dens),alpha=0.7)
	#points(ab[,'alpha'],ab[,'beta'],col=rgb(t(col2rgb(col)/255), alpha=.1),pch=pch,cex=3*cex)
	points(ab[,'alpha'],ab[,'beta'],col=col,pch=20,cex=cex)
}

