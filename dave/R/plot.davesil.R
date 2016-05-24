plot.davesil<- function(x,...,range=NULL) {
	o.davesil<- x
	test<- is.null(range) 
	ndim<- length(o.davesil$o.relgr)
	if(test == TRUE) range<- c(1,ndim)
	n<- range[2]-range[1]+1
# organize o.sil like output of function silhouette
	o.sil<- cbind(o.davesil$o.relgr,o.davesil$second,o.davesil$swidth)
#
	oo.sil<- o.sil[order(o.sil[,3]),]
	ooo.sil<-oo.sil[order(oo.sil[,1]),]
	ooo.sil<- ooo.sil[range[1]:range[2],]
	groups<- unique(ooo.sil[,1])
	ngroups<- length(groups)-1
	tgroups<- table(ooo.sil[,1])
	space<- rep(0,n)
	yc<- seq(1,n,1)
	pos<- 0
	if(ngroups > 0) {
		for (i in 1:ngroups) {
			pos<- pos+tgroups[i] 
			space[pos+1]<- 0.2
		}
	}
	for (i in 1:n) yc[i]<- yc[i]+sum(space[1:i])
#  cat("ngroups ",ngroups,"\n")
#  cat("groups ",groups,"\n")
#  cat("tgroups ",tgroups,"\n")
#  cat("space ",space,"\n")
#  cat("yc ",yc,"\n")
	n1<- o.davesil$names
	n2<- n1[order(o.sil[,3])]
	n3<- n2[order(oo.sil[,1])]
	n3<- n3[range[1]:range[2]]
	par(tcl=-0.3,mgp=c(2,0.5,0))
	barplot(ooo.sil[,3],space=space,horiz=TRUE,border=0,xlim=c(-0.4,1),xlab="Silhouette width",ylab="",col=gray(0.9/ooo.sil[,1]^0.5),cex.axis=0.8,cex.lab=0.8,xpd=TRUE,srt=0,axisnames=FALSE)
	legend("topright",table(ooo.sil[,1]),paste("group no.",unique(ooo.sil[,1])),pch=15,col=gray(0.9/(unique(ooo.sil[,1]))^0.5),bty="n",cex=0.8)
# These are the releve numbers on the left hand side of the barplot. Fontsize is a heuristic function of ndim:
	textsize<- 3.5/ndim^0.5
#	cat("textsize ",textsize,"\n")
	text(rep(-0.35,ndim),yc,n3,cex=textsize,pos=1,offset=0.0,adj=c(1,0))
}
