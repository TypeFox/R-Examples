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


pi<-get.pi(x,fun,r=r.max,r.low=r.min)
pi.null<-get.pi.permute(x,fun,r=r.max,r.low=r.min,permutations=100)

null.ci<-apply(pi.null,2,quantile,probs=c(0.25,0.75))

plot(r.mid, pi , type="l")
lines(r.mid, null.ci[1,] , lty=2)
lines(r.mid, null.ci[2,] , lty=2)
