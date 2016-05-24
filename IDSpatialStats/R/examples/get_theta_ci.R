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


theta<-get.theta(x,fun,r=r.max,r.low=r.min)
theta.ci<-get.theta.ci(x,fun,r=r.max,r.low=r.min,boot.iter=100)


plot(r.mid, theta , type="l")
lines(r.mid, theta.ci[1,] , lty=2)
lines(r.mid, theta.ci[2,] , lty=2)
