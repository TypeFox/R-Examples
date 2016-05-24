image.matrixBC <-
function(x,labelx=TRUE,labely=TRUE,...){
res <- x$mus
dev.new(width=6,height=6)
image(x=1:ncol(res),y=1:nrow(res),t(res),axes=TRUE,xaxt='n',yaxt='n',main="",xlab="",ylab="",col=heat.colors(12,alpha=0.8),...)

if(labelx==TRUE){
	axis(1,at=1:ncol(res),cex.axis=1.2)
}
if(labely==TRUE){
	axis(2,at=1:nrow(res),labels=c(nrow(res):1),cex.axis=1.2)
}

title(main="Estimated Mean Matrix",xlab="Features",ylab="Observations",cex.lab=1.3,cex.main=1.5)

invisible(x)
invisible(res)

}
