#to get a quick look at the consistent subsets at the 95% level after correction
#for multiple comparisons
#Plots the p-values for a z-test of |x[i]-x[j]|/sqrt(u[i]^2+u[j]^2) == 0, (two-tailed)
#key.height is a fraction of the figure region height
#key.width is the width of the key area in cm, unless under 1, 
#in which case it is interpreted as a fraction of the plot region width.

cplot<-function(x,u,labels=names(x), p.adjust.method="holm", ordered=TRUE,
	breaks=c(0,0.001,0.01, 0.05, 0.1,1), col=terrain.colors(length(breaks)-1), log.p=FALSE,
	main=paste("Consistency map -", deparse(substitute(x))), subtitle=NULL, key=FALSE, 
	key.width=2.54, key.height=0.6,...) {
	
	L<-length(x)
	ID<-1:L
	names(ID)<-as.character(labels)
	
	p.xx<-function(i,j,x,u) 
		2*pnorm(abs(x[i]-x[j]), 0, sqrt(u[i]^2+u[j]^2), lower.tail=FALSE)
	
	if(ordered) oo <- order(x) 
		else oo <- 1:L
	
	p<-outer(1:L, 1:L, FUN=p.xx,x=x[oo],u=u[oo])
	
	p[upper.tri(p)]<-p.adjust(p[upper.tri(p)], method=p.adjust.method)
	p[lower.tri(p)]<-p.adjust(p[lower.tri(p)], method=p.adjust.method)
	diag(p)<-rep(1, L) #Because it always is!
	
	rownames(p)<-colnames(p)<-labels[oo]
	
	if(log.p) p<- -log10(p)
	
	if(key) {
 		mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
		if(key.width<1) {
			layout(matrix(c(2, 1), ncol = 2), widths = c(1, key.width))
		} else {
			w <- (3 + mar.orig[2]) * par("csi") * key.width
			layout(matrix(c(2, 1), ncol = 2), widths = c(1, lcm(w)))
		}
		mar <- mar.orig
		mar[4] <- mar[2]
		mar[2] <- 1
		par(mar = mar)
 		par(mai=c((1-key.height)*par("fin")[2], par("mai")[2:4]))
		
		z<-matrix(breaks[-1]-diff(breaks)/2, nrow=1)
		
		image(z=z, breaks=breaks, col=col, xaxt="n", yaxt="n",xlab="", ylab="")
		yl<-par("usr")[3:4]
		yt<-seq(yl[1],yl[2], along=breaks)
		axis(4,at=yt, labels=breaks, las=1, cex.axis=0.7)
		box()
		if(log.p) mtext(expression(-log[10](p)), side=1,line=0.2)
			else mtext("p", side=1,line=0.2)
		
		par(mar=par.orig$mar)

	}

	
	image(x=1:L, y=1:L, z=p, breaks=breaks, col=col, xaxt="n", yaxt="n",xlab="", ylab="",...)	
	box()
	axis(1,at=1:L, labels=labels[oo], las=2)
	axis(2,at=1:L, labels=labels[oo], las=1)
	title(main=main)
	if(log.p && is.null(subtitle)) mtext(expression(-log[10](p)), side=3,line=0.2)
	
	return(invisible(p))
}


