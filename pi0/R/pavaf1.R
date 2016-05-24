pavaf1=function(p,max.bins=20,bin.method=c("max","Sturges","Scott","FD"),
                discrete=FALSE,plotit=FALSE,...)
{
    #library("Iso")
    if(discrete){
        counts=table(p)
        cents=as.numeric(names(counts))
        dens=counts/sum(counts)*length(counts)
        if(plotit)plot(cents,dens/length(counts),type='h',lwd=2,xlim=c(0,1),xlab='p',ylab='probability')
    }else{
        bin.method=match.arg(bin.method)
        if (bin.method=='max'){
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
            histobj=hist(p,breaks=brs,probability=TRUE,xlim=c(0,1),xlab='p',...)
            box()
        }else{
            histobj=hist(p,breaks=brs,plot=FALSE)
        }
        dens=histobj$density
        cents=histobj$mids
    }
    pavaed=pava(dens,decreasing=TRUE)
    if(plotit){
        lines(cents,pavaed/if(discrete) length(counts) else 1,col=4,lwd=2)
        abline(h=tail(pavaed,1)/if(discrete) length(counts) else 1,col=4,lwd=1)
    }
    tail(pavaed,1)
}
