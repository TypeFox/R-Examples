oneway.plots <-
function(Y,fac1,COL=c("#A9E2FF","#0080FF"))
{
par(mfrow=c(1,3),pty="s")
YL <- range(Y)
dotchart(x=Y,groups=fac1,col=COL[2],pch=1,
xlab=deparse(substitute(Y)),xlim=YL,
gdata=tapply(Y,fac1,mean),gpch=17)
plot(Y~fac1,col=COL[1],ylab=deparse(substitute(Y)),ylim=YL)
plot.design(Y~fac1,ylab=deparse(substitute(Y)),ylim=YL)
par(mfrow=c(1,1),pty="m")
}

