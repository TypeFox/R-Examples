plot.triang.list <- function(x,showrect=FALSE,add=FALSE,xlab="x",ylab="y",
                             main="",asp=1,...) {
stopifnot(inherits(x,"triang.list"))
rw <- attr(x,"rw")
if(!add) {
	plot(0,0,type="n",xlim=rw[1:2],ylim=rw[3:4],
             xlab=xlab,ylab=ylab,main=main,asp=asp)
}
for(tri in x) {
	polygon(as.list(tri),...)
}
if(showrect) do.call(rect,as.list(rw)[c(1,3,2,4)])
invisible()
}
