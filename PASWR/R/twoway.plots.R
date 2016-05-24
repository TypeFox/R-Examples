twoway.plots <-
function(Y, fac1, fac2, COL=c("#A9E2FF","#0080FF"))
{
 par(mfrow=c(2,2),pty="m",mar=c(5.1,4.1,4.1,2.1))
 YL <- range(Y)
 plot(Y~fac1,col=COL[1],main=deparse(substitute(fac1)),
 xlab="",ylab=deparse(substitute(Y)),ylim=YL)
 plot(Y~fac2,col=COL[2],main=deparse(substitute(fac2)),xlab="",
 ylab=deparse(substitute(Y)),ylim=YL)
 plot.design(Y~fac1+fac2,fun="mean",ylab=deparse(substitute(Y)),ylim=YL)
 interaction.plot(fac1,fac2,Y,xlab=deparse(substitute(fac1)),
 main="Interaction",trace.label=deparse(substitute(fac2)),type="b",
 legend=FALSE,ylab=deparse(substitute(Y)),ylim=YL)
 par(mfrow=c(1,1))
}

