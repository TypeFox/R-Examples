recluster.plot <-function (tree, data, low=1, high=0, id=NULL,nodelab.cex=0.8, direction="downwards",...) {

	if(length(tree$tip.label)!= (nrow(data)+1)) stop("ERROR: different site numbers between tree and bootstrap matrix")

	plot(tree, no.margin=T, direction=direction,...)
			

	if (direction=="downwards") {
		nodelabels(round(data[,low]),cex=nodelab.cex, adj = c(1.2, -0.3), frame="none", col=c(id))
		if(high>0) {
			nodelabels(round(data[,high]),cex=nodelab.cex, adj = c(-0.2, -0.3), frame="none", col=c(id))
		}
	}
	if (direction=="rightwards") {
		nodelabels(round(data[,low]),cex=nodelab.cex, adj = c(1.2, -0.3), frame="none", col=c(id))
		if(high>0) {
			nodelabels(round(data[,high]),cex=nodelab.cex, adj = c(1.2, +1.2), frame="none", col=c(id))
		}
	}
}
