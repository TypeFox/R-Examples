ks.sporm <-
function(x, y, B=1000){
    res<-ks.stat(x,y)
    ks<-res$ks
    theta<-res$theta
    m<-length(x)
    n<-length(y)
    Rej<-NULL
    Simul<-function(){
        u<-runif(m)
        v<-runif(n)
        v<-v/(theta-(theta-1)*v)
        res<-try(ks.stat(u, v), TRUE)
        (res$ks>ks)
    }
    Rej <- lapply(1:B, function(i) try(Simul(), TRUE))
    Rej<-unlist(Rej[sapply(Rej, function(x) !inherits(x, "try-error"))])
    pval<-mean(Rej)
    list(theta=theta,ks=ks, pval=pval)#, B=length(Rej))
}
