`agjack.pi0` <-
function(subtobj,mean.n=c('mean','harmean','geomean') ,pointpair=FALSE, trunc=TRUE,tol=1e-5){
    y=subtobj[,'f1']
    n1=subtobj[,'n1']
    n2=subtobj[,'n2']
    min.n=pmin(n1,n2)
    max.n=pmax(n1,n2)
    mean.n=match.arg(mean.n)
    mean.n=get(mean.n)
    n=apply(cbind(n1,n2),1,mean.n)
    
    if(!pointpair){
        y.mean=aggregate(y,list(min.n,max.n),mean)
        lev=apply(y.mean[,1:2],2,function(xx)as.numeric(xx))
        y.mean=y.mean[,3]
        n.unique=apply(lev,1,mean.n)
        nns=length(n.unique)
        rslt=0
        k=0
        for(i in 1:(nns-1))
            for(j in (i+1):nns)
                if(abs(sqrt((n.unique[j]+2)/(n.unique[i]+2))-1)>tol) {
                    k=k+1
                    rslt=rslt+if(trunc){
                        min(1,max(0,gjack(y.mean[i],y.mean[j],sqrt((n.unique[j]+2)/(n.unique[i]+2))),na.rm=TRUE),na.rm=TRUE)
                        }else{gjack(y.mean[i],y.mean[j],sqrt((n.unique[j]+2)/(n.unique[i]+2)))}
                }
        return(rslt/k)
    }else{
        nns=length(n)
        rslt=0
        k=0
        for(i in 1:nns)
            for(j in 1:nns)
                if(abs(sqrt((n[j]+2)/(n[i]+2))-1)>tol) {
                    k=k+1
                    rslt=rslt+if(trunc){
                            min(1,max(0,gjack(y[i],y[j],sqrt((n[j]+2)/(n[i]+2))),na.rm=TRUE),na.rm=TRUE)
                            }else{gjack(y[i],y[j],sqrt((n[j]+2)/(n[i]+2)))}
                }
        return(rslt/k)
    }
}

