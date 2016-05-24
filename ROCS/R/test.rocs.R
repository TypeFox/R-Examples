test.rocs <-
function(x0, x1, z0, z1, B=1000, do.plot=TRUE)
{
    vus.x<-rocs.x(x0,x1,n.perm=1,do.plot=FALSE)[[1]]
    vus.z<-rocs.x(z0,z1,n.perm=1,do.plot=FALSE)[[1]]
    vus.d<-vus.x-vus.z
    
    l.x0<-length(x0)
    l.x1<-length(x1)
    l.z0<-length(z0)
    l.z1<-length(z1)
    
    target.l.0<-round((l.x0+l.z0)/2)
    target.l.1<-round((l.x1+l.z1)/2)
    
    my.oversample<-function(a,n)
    {
        la<-length(a)
        if(la>=n)
        {
            b<-sample(a, n, replace=FALSE)
        }else{
            reps<-floor(n/la)
            b<-c(rep(a, reps), sample(a, n-la*reps, replace=FALSE))
        }
        b
    }
    
    x0<-my.oversample(x0, target.l.0)
    x1<-my.oversample(x1, target.l.1)
    z0<-my.oversample(z0, target.l.0)
    z1<-my.oversample(z1, target.l.1)
    
    new.x<-c(x0, x1)
    new.y<-c(rep(0, length(x0)), rep(1, length(x1)))
    new.x<-(rank(new.x, ties.method = "random")-1)/(length(x0)+length(x1)-1)
    new.x0<-new.x[new.y == 0]
    new.x1<-new.x[new.y == 1]
    
    new.z<-c(z0, z1)
    new.y<-c(rep(0, length(z0)), rep(1, length(z1)))
    new.z<-(rank(new.z, ties.method = "random")-1)/(length(z0)+length(z1)-1)
    new.z0<-new.z[new.y == 0]
    new.z1<-new.z[new.y == 1]
    
    #
    
    merged.0<-c(new.x0, new.z0)
    merged.1<-c(new.x1, new.z1)
    
    ref.distr<-1:B
    for(i in 1:B)
    {
        this.x0<-sample(merged.0, l.x0, replace=TRUE)
        this.x1<-sample(merged.1, l.x1, replace=TRUE)
        this.z0<-sample(merged.0, l.z0, replace=TRUE)
        this.z1<-sample(merged.1, l.z1, replace=TRUE)
        
        this.vus.x<-rocs.x(this.x0, this.x1,n.perm=1,do.plot=FALSE)[[1]]
        this.vus.z<-rocs.x(this.z0, this.z1,n.perm=1,do.plot=FALSE)[[1]]
        
        ref.distr[i]<-this.vus.x - this.vus.z
    }
    pval<-min(sum(ref.distr >= vus.d)/B, sum(ref.distr <= vus.d)/B)
    
    if(do.plot)
    {
        r<-range(c(ref.distr, vus.d))
        hist(ref.distr, xlim=r)
        abline(v=vus.d, col="red",lwd=2)
    }
    
    return(2*pval)
}
