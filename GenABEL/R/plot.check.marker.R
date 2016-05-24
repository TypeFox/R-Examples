"plot.check.marker" <- 
function(x, y, ...) {
	if (!is(x,"check.marker")) stop("wrong class of x (should be produced by check.marker()")
	xx <- c()
	y <- c()
	g <- c()
	gid <- 1
	lims <- c(min(x$call$map),max(x$call$map))
	plot(lims[1],lims[1],col=NA,xlim=lims,ylim=lims,xlab="",ylab="", ...)
    if (length(x$details.redundancy)>1) {
	for (i in 1:length(x$details)) {
		namx <- names(x$details.redundancy[i])
		if (namx == "all") next;
		namy <- x$details[[i]]
		cooy <- x$call$map[match(namx,x$call$name)]
		text(x=cooy,y=cooy,labels=c(namx),adj=c(1,0),col="#000000",cex=.75)
		xx <- c(xx,cooy)
		exty <- x$call$map[match(namy,x$call$name)]
		if (length(exty) == 0) next
		xx <- c(xx,exty)
		y <- c(y,rep(cooy,length(exty)+1))
		g <- c(g,rep(gid,length(exty)+1))
		gid <- gid + 1
	}
	colr <- c("#000000","#FF0000","#00FF00","#0000FF","#FF00FF","#00FFFF","#FFFF00","#990000","#009900","#000099","#990099","#009999","#999900")
	if (length(g) >0) for (i in 1:max(g)) {
		j <- i %% length(colr)
		if (j == 0) j = length(colr)
		col <- colr[j]
#		abline(h=y[g==i][1],col=col,lty=2)
		lines(c(xx[g==i][1],xx[g==i][length(xx[g==i])]),c(xx[g==i][1],xx[g==i][1]),col=col,lty=1)
		lines(c(xx[g==i][1],xx[g==i][1]),c(xx[g==i][1],0),col=col,lty=1)
		for (j in 2:length(xx[g==i])) {
#			abline(v=xx[g==i][j],col=col,lty=2)
			lines(c(xx[g==i][j],xx[g==i][j]),c(xx[g==i][1],0),col=col,lty=2)
		}
		points(xx[g==i][2:length(xx[g==i])],y[g==i][2:length(xx[g==i])],col=col,pch=20,cex=1.)
		points(xx[g==i][1],y[g==i][1],col=col,pch=3,cex=1.)
	}
    }
	abline(0,1)
	rug(x$call$map)
	cred <- match(x$redundant,x$call$name) ###
	cred <- cred[!is.na(cred)]
	if (length(cred)>0) {
		coo <- x$call$map[cred]
		rug(coo,col="cyan")
	}
	cred <- match(x$nohwe,x$call$name) ###
	cred <- cred[!is.na(cred)]
	if (length(cred)>0) {
		coo <- x$call$map[cred]
		rug(coo,col="green")
	}
	cred <- match(x$nofreq,x$call$name) ###
	cred <- cred[!is.na(cred)]
	if (length(cred)>0) {
		coo <- x$call$map[cred]
		rug(coo,col="yellow")
	}
	cred <- match(x$nocall,x$call$name) ###
	cred <- cred[!is.na(cred)]
	if (length(cred)>0) {
		coo <- x$call$map[cred]
		rug(coo,col="red")
	}
	rug(c(min(x$call$map),max(x$call$map)),col="black")
	cat("Red: no call\nYellow: low frequency\nGreen: out of HWE\n")
	cat("Cyan: redundant\nDiagonal: redundant markers (reference presented as \"+\")\n")
}

