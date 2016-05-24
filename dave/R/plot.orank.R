plot.orank<- function(x,...) {
	plot(x$x.axis,x$y.axis,asp=1,xlab="",ylab="",axes=FALSE,type="n")
	points(x$x.axis,x$y.axis,pch=16,col="black",cex=0.9)
	if(x$use == "columns") allnam<- substr(x$all.colnam,1,20)
	if(x$use == "rows") allnam<- x$all.rownam
	text(x$x.axis,x$y.axis,allnam,cex=0.4,col="black",pos=4)
# mark the 5 plots ranked highest
	for(i in 1:x$n.ranks) {
		j<- which(allnam == x$var.names[i])
		points(x$x.axis[j],x$y.axis[j],pch=16,col="black",cex=1.6)
		text(x$x.axis[j],x$y.axis[j],as.character(i),cex=1.0,pos=3,offset=0.4)
	}
}
