plot.pdclust <-
function(x,labels=NULL,type="rectangle",cols="black", timeseries.as.labels = T, p.values=F, ...)
{
	# create labels if none are given
	#if (is.null(labels)) {
#		labels <- paste("",1:length(X),sep="");
#	}

	if (!is.null(labels)) {
		timeseries.as.labels=F
	}
	
	# create color vector for plotting ts
	if (length(cols) == 1) {
		cols = rep(cols, x$N);
	}

	
	if (!timeseries.as.labels )
	{
		class(x) <- "hclust"
		plot(x, main="Permutation Distribution Clustering",labels, ...);
		
		if (p.values) {
			plot.add.pvalues(x)
		}
	} else {
		
		X <- x$data

			# -- or horizontal plot --
		oldpar <- par(no.readonly=TRUE)
		par(fig=c(0,0.3,0,1),mar=c(0,0,0,0),new=F)
		on.exit(par(oldpar))
		class(x) <- "hclust"
		plot(as.dendrogram(x), 
		horiz = TRUE, leaflab="none", frame.plot=FALSE,
		 edgePar=list(col='black',lw=2),
		 type=type,...
		 )
	

		scale <- 0.9
		offset <- 0.05
		for (i in 1:x$N) {

			par(fig=c(0.35,1,offset+(i-1)/x$N*scale,offset+(i/x$N*scale)),mar=c(0,0,0,0),new=TRUE)
			if (x$multichannel) {
				dat <- X[,x$order[i],1]
			} else {
				dat <- X[,x$order[i]]
			}
			 plot(dat,ann=F,xaxt='n',yaxt="n",frame.plot=F,type="n")
			 lines(dat,lw=2,col=cols[x$order[i] ])

			#text( 0, 0,  labels[i])
			
			
		}
		
		
	}
	
			
	invisible()
	
}

sw<-function(pic)
{
	temp <- (pic[,,1]+pic[,,2]+pic[,,3])/3
	pic[,,1] <- temp
	pic[,,2] <- temp
	pic[,,3] <- temp
	return(pic)
}

