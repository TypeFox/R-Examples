betagen=function(Nsim=10^3){
# Generation from a Beta(a,b) target
# using a uniform and a Beta([a],[b]) proposals
a=2.7;b=6.3

M=2.67
y=runif(Nsim)
u=runif(Nsim,max=M)
par(mfrow=c(1,2),mar=c(4,4,2,1))
plot(y,u,col="grey",pch=19,cex=.4,ylab=expression(u.g(y)))
points(y[u<dbeta(y,a,b)],u[u<dbeta(y,a,b)],pch=19,cex=.4)
curve(dbeta(x,a,b),col="sienna",lwd=2,add=T)
abline(h=M,col="gold4",lwd=2)

M=1.68
y=rbeta(Nsim,2,6)
u=runif(Nsim,max=M)
labels=u<dbeta(y,a,b)/dbeta(y,2,6)
plot(y,u*dbeta(y,2,6),col="grey",pch=19,cex=.4,ylab=expression(u.g(y)))
points(y[labels],u[labels]*dbeta(y[labels],2,6),pch=19,cex=.4)
curve(dbeta(x,a,b),col="sienna",lwd=2,add=T)
curve(M*dbeta(x,2,6),col="gold4",lwd=2,add=T)
}
