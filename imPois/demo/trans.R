# 
# This demo illustrates a translation behaviour of imprecise prior 
# when a sample is newly obsered.  This geometrical representation
# gives a convenient approach of finding maximum and minimum. 
# 
# updated on 2015.11.27

e1 <- seq(from=0.1, to=10, by=0.2) # shape 
e0 <- seq(from=0.1, to=10, by=0.2) # rate 
nps <- expand.grid(e1=e1,e0=e0)
nps$et <- apply(nps, 1, function(x) evfn(y=numeric(0), pars=c(0,x[1],x[2]))$value)


# 2d: log-gamma (ztrunc=FALSE)
lambda <- 1
n <- 1e1
set.seed(16979238) # 1,1,1,0,0,2,1,0,2,
y <- rpois(n=n, lambda=lambda)

contour(x=e1, y=e0, z=matrix(nps$et, nrow=length(e1), ncol=length(e0)), levels=seq(from=-3, to=3, by=0.5), method="edge", labcex=1.5, lty="dashed", xlab=expression(xi[1]), ylab=expression(xi[0]), cex.lab=1.5, xlim=c(0,10), ylim=c(0,10))
polygon(x=c(0,0,1,1),y=c(0,1,1,0), border="darkblue")
polygon(x=c(0,1),y=c(1,0), lwd=2)

convh <- iprior(ui=rbind(diag(2), -diag(2)), ci=c(0,0,-1,-1)) 
for(i in 1:n){
	op <- update(convh, y=y[1:i], ztrunc=FALSE)
	polygon(op$vtx1)
	vtx1.m <- colMeans(op$vtx1)
	text(x=0.2+vtx1.m[1], y=0.2+vtx1.m[2], labels=bquote(y[.(i)]), cex=1.5)
	sop <- summary(op)
	text(x=(sop$inf.p1)[1], y=(sop$inf.p1)[2], labels=round(sop$inf,3), col="darkblue", cex=0.75)
	text(x=(sop$sup.p1)[1], y=(sop$sup.p1)[2], labels=round(sop$sup,3), col="red", cex=0.75)
}
legend("bottomright", legend=c("min", "max"), col=c("darkblue", "red"), lty=c(1,1), lwd=c(2,2))

