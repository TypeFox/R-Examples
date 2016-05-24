#class: usually disease labels
#class2: different sources like multiple studies
PlotPC2D <- function(coord, drawObjects=TRUE, drawEllipse=FALSE, dataset.name=NULL, 
		pctInfo=NULL, main=NULL, sub=NULL, xlab=NULL, ylab=NULL, newPlot=TRUE, .points.size=1,
		.class=rownames(coord), .class.order=NULL, .class.color=NULL, .class2=NULL, .class2.order=NULL, .class2.shape=NULL, 
		.annotation=TRUE, .legend=c("bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right", "center")) 
{
	.legend <- match.arg(.legend)
	if(is.null(.class.order)) {
		.label <- factor(.class) 		
	} else {
		.label <- factor(.class, levels=.class.order)
	}
	if(is.null(.class.color)) {
		.class.color <- 1:nlevels(.label) 		
	} 
	
	if(!is.null(.class2)) {
		if(is.null(.class2.order)) {
			.label2 <- factor(.class2) 		
		} else {
			.label2 <- factor(.class2, levels=.class2.order)
		}
		if(is.null(.class2.shape)) {
			.class2.shape <- 1:nlevels(.label2)  		
		} 
	}
	
	colnames(coord) <- c("PC1","PC2")
	
	if(newPlot) {
		.sysname <- Sys.info()["sysname"]
		if(.sysname=="Windows")
			windows()
		else if(.sysname=="Darwin")
			quartz()
		else #*nix
			try(X11(), TRUE)
	}
	
	
	if(drawEllipse) {
		requireAll(c("ellipse","MASS"))
		.ellip <- NULL 
		for(lev in levels(.label)) {
			.mve <- cov.mve(coord[.label==lev,1:2], method="classical")
			.ellip <- rbind(.ellip, ellipse(x=.mve$cov, centre=.mve$center, npoints=100))
		}
		plot(rbind(.ellip,coord), type="n", xlab="", ylab="", axes=FALSE, xlim=range(rbind(.ellip,coord)), ylim=range(rbind(.ellip,coord)))		
		for(i in 1:nlevels(.label)) {
			lines(.ellip[(1+100*(i-1)):(100*i),], col=.class.color[i])
		}
	} else {
		plot(coord, type="n", xlab="", ylab="", axes=FALSE, xlim=range(coord), ylim=range(coord))	
	}
	
	if(drawObjects) 
		for(i in 1:nlevels(.label)) {
			if(exists('.label2')) {
				for(j in 1:nlevels(.label2)) 
					points(coord[.label==levels(.label)[i] & .label2==levels(.label2)[j],], pch=.class2.shape[j], col=.class.color[i], cex=.points.size)
			} else
				points(coord[.label==levels(.label)[i],], pch=20, col=.class.color[i], cex=.points.size)			
		}
	
	if(.annotation) {
		axis(side=1,lwd=4,tck=-0.02)
		axis(side=2,lwd=4,tck=-0.02)
	}
	box(bty="L",lwd=4)
	abline(v=0, lty=2)
	abline(h=0, lty=2)	
	
	if(.annotation) {
		if(exists('.label2')) {
			legend2(.legend, legend=c(levels(.label),levels(.label2)), pch=c(rep(-1,nlevels(.label)), .class2.shape), fill=c(.class.color, rep(-1,nlevels(.label2))))
		} else
			legend(.legend, legend=levels(.label), text.col=.class.color, pch=20, col=.class.color)
		
		.pctInfo <- c(ifelse(is.null(pctInfo), "", paste(" (",pctInfo[1],"%)",sep="")),
				ifelse(is.null(pctInfo), "", paste(" (",pctInfo[2],"%)",sep="")))
		
		.main <- ifelse(is.null(main), paste("Meta-PCA [", dataset.name, "]", sep=""), main)
		.xlab <- ifelse(is.null(xlab), paste("Meta-PC1", .pctInfo[1], sep=""), xlab)
		.ylab <- ifelse(is.null(ylab), paste("Meta-PC2", .pctInfo[2], sep=""), ylab)
		
		title(main=.main, sub=sub, xlab=.xlab, ylab=.ylab)
	}
	invisible()
}
