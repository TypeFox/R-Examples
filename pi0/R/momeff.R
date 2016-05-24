parncpt.momeff=function(tstat,n1,n2=n1,zeromean,gamma2,lower.df=6.1,upper.df=100,approx=TRUE) {
###### did not deal with missing n1 case yet!!!

df=v=n1+n2-2

if(missing(zeromean) ) {
    if(n1+n2-2<=4)zeromean=FALSE else zeromean='auto'
}

if(n1+n2-2<=4 )zeromean=FALSE

if(zeromean=='auto'){

    ans.zeromean=if(missing(gamma2)){
            momeff.zeromean(tstat,n1,n2,,lower.df,upper.df,approx)
    }else   momeff.zeromean(tstat,n1,n2,gamma2,lower.df,upper.df,approx)

    ans.nonzeromean=if(missing(gamma2)){
            momeff.nonzeromean(tstat,n1,n2,,lower.df,upper.df,approx)
    }else   momeff.nonzeromean(tstat,n1,n2,gamma2,lower.df,upper.df,approx)

    if(attr(ans.zeromean,'gamma2') * attr(ans.nonzeromean,'gamma2')>0){
        ## either both have positive gamma2, or neither has
        ans=c(ans.zeromean,ans.nonzeromean)
        ncp.mean=c(0,attr(ans.nonzeromean,'eff.mean')/sqrt(1/n1+1/n2))
        ncp.var=c(attr(ans.zeromean,'gamma2')/(1/n1+1/n2),attr(ans.nonzeromean,'gamma2')/(1/n1+1/n2))
        if(any(ncp.var<0))ncp.var=c(1e-3,1e-3)

        quants=quantile(tstat,seq(0,1,length=500))
        deltas=seq(quantile(tstat,.001),quantile(tstat,.999),length=201)
        dd=suppressWarnings(outer(quants,(deltas[-1]+deltas[-length(deltas)])/2,
                function(xx,yy)dt(xx,v,yy)))
        ll=c(-sum(log(
                ans[1]*dt(quants,v)+(1-ans[1])*dd%*%(pnorm(deltas[-1],ncp.mean[1],sqrt(ncp.var[1]))-
                                        pnorm(deltas[-length(deltas)],ncp.mean[1],sqrt(ncp.var[1])))
                )),-sum(log(
                ans[2]*dt(quants,v)+(1-ans[2])*dd%*%(pnorm(deltas[-1],ncp.mean[2],sqrt(ncp.var[2]))-
                                        pnorm(deltas[-length(deltas)],ncp.mean[2],sqrt(ncp.var[2])))
                )))
        ans=if(which.min(ll)==1) ans.zeromean else ans.nonzeromean
    }else {## one has positive gamma2, the other does not
        ans=if(attr(ans.zeromean,'gamma2')>0) ans.zeromean else ans.nonzeromean
    }

}else if (zeromean){
    ans.zeromean=if(missing(gamma2)){
            momeff.zeromean(tstat,n1,n2,,lower.df,upper.df,approx)
    }else   momeff.zeromean(tstat,n1,n2,gamma2,lower.df,upper.df,approx)
    ans=ans.zeromean
}else {
    ans.nonzeromean=if(missing(gamma2)){
            momeff.nonzeromean(tstat,n1,n2,,lower.df,upper.df,approx)
    }else   momeff.nonzeromean(tstat,n1,n2,gamma2,lower.df,upper.df,approx)
    ans=ans.nonzeromean
}

ans
}



