histf1=function(p,max.bins=20,bin.method=c("max","nmse","bootstrap","Sturges","Scott","FD"),
        discrete=FALSE,seq.perm=FALSE,nboots=200,rightBoundary=FALSE, plotit=FALSE,perm.n,perm.h,...)
{
    if(seq.perm){
        if(missing(perm.h) || missing(perm.n)){
            vals=c(0,sort(unique(p)))
        }else{
            vals=c(0:(perm.h-1)/perm.n,perm.h/(perm.n:perm.h))
        }
        p=factor(p,levels=vals[-1])
        wids=diff(vals)
        cts=table(p)
        dens=cts/sum(cts)*length(vals[-1])
        if(plotit)plot(vals[-1],dens,type='h',lwd=2,xlim=c(0,1),xlab='p',ylab='density')
    }else if(discrete){ 
        if(missing(perm.n)){
            vals=c(0,sort(unique(p)))
        }else{
            vals=0:perm.n/perm.n
        }
        p=factor(p,levels=vals[-1])
        cts=table(p)
        dens=cts/sum(cts)*length(vals[-1])
        wids=rep(1,length(vals[-1]))/length(vals[-1])
        if(plotit)plot(vals[-1],dens,type='h',lwd=2,xlim=c(0,1),xlab='p',ylab='density')
    }else{
        bin.method=match.arg(bin.method)
        ff=function(bb,pp){
            dens=hist(pp,br=seq(0,1,length=bb+1),plot=FALSE)$dens
            tailmean=rev(cumsum(rev(dens))/seq_along(dens))
            tailmean[which.max(dens<=tailmean)]
        }
        if(bin.method=='bootstrap'){
            min.pi0=min(sapply(2:max.bins,ff,pp=p))
            bootpi0s=matrix(,nboots,max.bins-1)
            f=function(bb,...)histf1(max.bins=bb,...)
            for(r in 1:nboots){
                pp=sample(p,replace=TRUE)
                bootpi0s[r,]=sapply(2:max.bins,ff,pp=pp)
            }
            n.bins=which.min(colMeans((bootpi0s-min.pi0)^2))+1
            brs=seq(0,1,length=n.bins+1)
        }else if (bin.method=='nmse'){
            pi0s=sapply(2:max.bins,ff,pp=p)
            nmse=(pi0s-min(pi0s))^2+pi0s*(1-pi0s)/length(p)
            n.bins=which.min(nmse)+1
            brs=seq(0,1,length=n.bins+1)
        }else if (bin.method=='max'){
            n.bins=max.bins
            brs=seq(0,1,length=n.bins+1)
        }else if (bin.method=='Sturges'){
            x=qnorm(p)
            n.bins=nclass.Sturges(x[is.finite(x)])
            rgx=range(x,na.rm=TRUE)
            brs=pnorm(seq(rgx[1],rgx[2],length=n.bins+1))
            brs=c(0,brs[-c(1,n.bins+1)],1)
        }else if (bin.method=='Scott'){
            x=qnorm(p)
            n.bins=nclass.scott(x[is.finite(x)])
            rgx=range(x,na.rm=TRUE)
            brs=pnorm(seq(rgx[1],rgx[2],length=n.bins+1))
            brs=c(0,brs[-c(1,n.bins+1)],1)
        }else if (bin.method=='FD'){
            x=qnorm(p)
            n.bins=nclass.FD(x[is.finite(x)])
            rgx=range(x,na.rm=TRUE)
            brs=pnorm(seq(rgx[1],rgx[2],length=n.bins+1))
            brs=c(0,brs[-c(1,n.bins+1)],1)
        }

        if(plotit){
            histobj=hist(p,breaks=brs,prob=TRUE,xlab='p',xlim=c(0,1),...)
        }else{
            histobj=hist(p,breaks=brs,plot=FALSE)
        }
        cts=histobj$counts
        wids=diff(histobj$breaks)
        dens=histobj$density
    }
    tailmean=rev(cumsum(rev(cts))/cumsum(rev(wids)))/length(p)
    f1=tailmean[min(which.max(c(dens<=tailmean, TRUE))+rightBoundary, length(cts))]
    if(plotit)abline(h=f1)
    drop(f1)
}

