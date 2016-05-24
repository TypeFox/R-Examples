trendplot <-
function(f,g, cex = 1,  cex.lab = 0.8, main=" ",prop=0.5,posleg="topleft",xlab="First Principal Axis",ylab="Second Principal Axis")
{
#------------------------------------------------------------------------------
# f and g are the row and column coordinates
#------------------------------------------------------------------------------
nrows<-dim(g)[[1]]
ncols<-dim(g)[[2]]
leg.txt<-dimnames(g)[[1]]
colsymb<-c(1:nrows)
gt<-t(g)
plot(f,g[1,],type="b", ylim = range(gt[1:ncols,],g[1:nrows,])/prop, xlab = xlab, ylab = ylab, cex = cex, cex.lab = cex.lab, main=main, col=1)
abline(h=0,lty=3)
   for (i in 1:(nrows)){
lines(f,g[i,],type="b",pch=i,  col=i)
}
legend(x=posleg,legend=leg.txt,col=colsymb,pch=c(1:nrows),bty="o",cex=.8)

}
