## This demo considers quantile regression spline estimation using
## crs() for a heteroskedastic Gaussian DGP.

library(crs)

set.seed(24)

n <- 10000

x <- runif(n,-1,1)
x.seq <- seq(min(x),max(x),length=100)

sd.x <- .1+abs(x)
sd.x.seq <- .1+abs(x.seq)

y <- rnorm(n,mean=x,sd=sd.x)

taus <- c(0.05,0.1,.25,0.5,.75,.9,0.95)


j.opt <- numeric()
s.opt <- numeric()

plot(x,y,cex=.25,col="lightgrey")
for(t in seq(along=taus)) {
  model.rq <- crs(y~x,tau=taus[t],cv.threshold=0,nmulti=2)
  lines(x.seq,predict(model.rq,newdata=data.frame(x=x.seq)),col=t,lty=t,lwd=2)
  lines(x.seq,qnorm(taus[t],mean=x.seq,sd=sd.x.seq),lty=2)
  j.opt[t] <- model.rq$degree
  s.opt[t] <- model.rq$segments
}

legend(min(x),max(y),paste("tau=",taus,", d=",j.opt,", s=",s.opt,sep=""),
       lty=1:length(taus),
       col=1:length(taus),
       lwd=rep(2,length(taus)),
       bty="n")



