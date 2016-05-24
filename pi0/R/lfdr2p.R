if(FALSE){

obj=function(par)-sum(log(abs(
    ((y-par[2])/(1-par[2])*beta(1,par[1]))^(1/(par[1]-1))
    /
    (par[1]-1)/(y-par[2])
    )))

lfdr2p=function(lfdr,pi0)
{
    y=pi0/lfdr
    fp1=min(y)
    phat=numeric(length(y))
        f.logy.fit=density(log(y),from=log(fp1),to=log(max(y)),n=length(y))
        f.logy.fun=approxfun(f.logy.fit$x,f.logy.fit$y)
    obj=function(b)abs(1/(y[i]-fp1)^((b-2)/(b-1))*(beta(1,b)/(1-fp1))^(1/(b-1))/(b-1)-fy)
    for(i in 1:length(y)){
        fy=f.logy.fun(log(y[i]))/y[i]
        b=nlminb(2,obj,lower=1+1e-2)$par
        phat[i]=1-((fp1-y[i])/(fp1-1)*beta(1,b))^(1/(b-1))
    }
    phat
}
#debug(lfdr2p)
z=qnorm(p/2)*c(-1,1)
tmp=locfdr(z,nulltype=0)
#tmp=locfdr(qnorm(p/2)*sample(c(-1,1),length(p),repl=T),nulltype=2)
y0=1/tmp$fdr*tmp$fp0[1,3]
phat=lfdr2p(tmp$fdr,tmp$fp0[1,3])


lfdr2p=function(lfdr,pi0)
{
    y=pi0/lfdr
    ord.y=order(y)
    G=length(y)
    fp1=min(y)
#    y[ord.y][1:sum(y==fp1)]=seq(fp1,min(setdiff(y,fp1)),length=sum(y==fp1))
    unif.idx=y==fp1
    n.fp1=sum(unif.idx)
    p.c=1-n.fp1/fp1/G
    sort.y=sort(y)
#    f.logy.fit=density(log(y),bw='ucv',kernel='epanechnikov',from=log(fp1),to=log(max(y)),n=length(y))
#    f.logy.fun=approxfun(f.logy.fit$x,f.logy.fit$y)
#    finv.der=f.logy.fun(log(y[ord.y]))/y[ord.y]^2
#    finv.der=dlogspline(log(y[ord.y]),
#            logspline(log(y),log(fp1*.99),log(max(y)*1.01)))/y[ord.y]^2
    dummy.y=seq(fp1,max(y),length=G*10)
    finv.der=dlogspline(dummy.y,
            logspline(y,(fp1*.99),(max(y)*1.01)))/dummy.y
    ints=(finv.der[-1]+finv.der[-length(finv.der)])*diff(dummy.y)/2
    ans=c(p.c,p.c-cumsum(ints/sum(ints)*p.c))
    ans=approx(dummy.y,ans,y[ord.y],rule=2)$y
#    finv.der=pmden(c(fp1,sort.y))$y[(n.fp1+1):G]
#    ints=(finv.der[-1]+finv.der[-length(finv.der)])*diff(y[ord.y[!unif.idx]])/2
#
#
#    ans=c(p.c,p.c-cumsum(ints))               ## theroretical
#    ans=c(p.c,p.c-cumsum(ints/sum(ints)*p.c))  ## intuitive surgery
    phat=numeric(G)
    phat[ord.y]=ans
#    phat[ord.y[!unif.idx]]=ans
    fact=function(fac,reps=100)
                mean(replicate(reps,{
                    phat[unif.idx]=runif(n.fp1,fac*p.c,1);
                    sum(pmax(0,diff(hist(phat,br=20,plot=FALSE)$density)))
                    }))
    decrease=sapply(31:50/50,fact)
    fac=(which.min(decrease)+30)/50
    phat[unif.idx]=runif(n.fp1,fac*p.c,1)
    phat
}

}
