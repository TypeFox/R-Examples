plot.FIT <-
function(fs,
	 plot.xlab="Gene Index",plot.ylab="FIT",plot.main="Main Title",
	 qqplot.ylab="FIT^2",qqplot.xlab="Chisq Quantitle",
	 pch.nonsig=21,pch.sig=19,col.pos="red",col.neg="blue")
    ## fs from feature.selection
{
    sx <- sort(qchisq(ppoints(fs$fullList$rd),df=1))
    sy <- sort(fs$fullList$rd)
    sy.sign <- fs$fullList$positive[order(fs$fullList$rd)]
    sy.significances <- fs$fullList$genes[order(fs$fullList$rd)]%in%fs$selectedList$genes
    plot(sx,sy,xlab=qqplot.xlab,ylab=qqplot.ylab,type="n")
    points(sx[!sy.significances],sy[!sy.significances],pch=pch.nonsig)
    points(sx[sy.significances&sy.sign],
           sy[sy.significances&sy.sign],
           pch=pch.sig,
           col=col.pos)
    points(sx[sy.significances&!sy.sign],
           sy[sy.significances&!sy.sign],
           pch=pch.sig,
           col=col.neg)
    abline(0,1,col="gray25",lty=2)
}
