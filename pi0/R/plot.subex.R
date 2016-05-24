`plot.subex` <-
function(x,y,rgl=TRUE,...)
{
    if(missing(y)) y=rgl else rgl=y
    if(rgl)rgl=requireNamespace("rgl",quitely=TRUE)
    if(!rgl) {dev.new(width=14,height=7);par(mfrow=c(1,2))}
    op=par( mar=c(5.1,4.1,4.1,4.1))
    hs=hist(x$pvalues,probability=TRUE,breaks=30,xlab='p-value',ylab='',
           ,border=4,axes=FALSE,main='Subsampling-Extrapolation based\np-value and q-value summary')
    box()
    height=max(hs$density)
    lines(c(-1,1),rep(x$pi0,2),col=4,lwd=2)
    mtext(expression(hat(pi)[0]),side=2,at=x$pi0,col=4,las=2,adj=2)


    ord=order(x$pvalues)
    lines(c(x$pvalues[ord],1,2),c(x$qvalues[ord],rep(x$pi0,2))*height,col=2,lwd=2)
    axis(1,at=seq(0,1,length=21),labels=rep('',21))
    axis(1,at=seq(0,1,length=11),labels=0:10/10)
    mtext(expression(hat(pi)[0]),side=4,at=x$pi0*height,col=2,las=2,adj=-3)

    axis(2,at=round(seq(0,height,length=8),1),col=4,col.axis=4)
    title(ylab='density',col.lab=4)

    axis(4,at=seq(0,height,length=21),labels=rep('',21),col=2,col.axis=2)
    axis(4,at=seq(0,height,length=11),labels=seq(0,1,length=11),col=2,col.axis=2)
    mtext('q-value',side=4,padj=3,col=2)

    par(op)
    plot(x$extrp.fit,rgl=rgl,...=...)
    invisible(NULL)
}

