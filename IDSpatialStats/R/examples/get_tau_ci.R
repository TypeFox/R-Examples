#compare normally distributed with uniform points
x<-cbind(1,runif(100,-100,100), runif(100,-100,100))
x<-rbind(x, cbind(2,rnorm(100,0,20), rnorm(100,0,20)))
colnames(x) <- c("type","x","y")

fun<-function(a,b) {
    if(a[1]!=2) return(3)
    if (b[1]==2) return(1)
    return(2)
}

r.max<-seq(10,100,10)
r.min<-seq(0,90,10)
r.mid <- (r.max+r.min)/2

tau<-get.tau(x,fun,r=r.max,r.low=r.min)
tau.ci<-get.tau.ci(x,fun,r=r.max,r.low=r.min,boot.iter=50)


plot(r.mid, tau ,ylim=c(1/max(tau.ci),max(tau.ci)), type="l", log="y")
lines(r.mid, tau.ci[1,] , lty=2)
lines(r.mid, tau.ci[2,] , lty=2)
lines(c(0,100),c(1,1), lty=3, col="grey")
