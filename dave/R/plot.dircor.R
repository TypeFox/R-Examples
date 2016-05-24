plot.dircor <-
function(x,...){
    o.dircor<- x
    step<- o.dircor$steps
    rlow<- o.dircor$lower.limit
    rupp<- o.dircor$upper.limit
    rout<- o.dircor$mean.correlation
#
    plot(c(step,step),c(rlow,rupp),xlab="Direction (deg)",ylab="Mantel correlation",type="n")
    points(step,rout,pch=1,cex=1.0)
    lines(step,rout,lwd=0.5)
    lines(step,rlow,lwd=0.3)
    lines(step,rupp,lwd=0.3)
    abline(v=c(50,100,150),lwd=1,col="gray")
    abline(h=0,lwd=1,col="gray")
 }
