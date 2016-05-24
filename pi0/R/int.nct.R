mTruncNorm=function(r=1, mu=0, sd=1, lower=-Inf, upper=Inf, approximation=c('int2','laplace', 'numerical'),integral.only=FALSE,...)
{   ## return E(X^r|lower<X<upper), where X~Normal(mu, sd), integral.only=FALSE (default); 
    ## return int_lower^upper X^r exp(-0.5*(X-mu)^2/sd^2) dx, if integral.only=TRUE
    if(any(lower>upper)) stop("boundary inconsistent")
    approximation=match.arg(approximation)
    if(integral.only) fact=1 else fact=1 /  (pnorm(upper,mu,sd)-pnorm(lower,mu,sd))/sqrt(2*pi)/sd 

    fact *
    if(approximation=='laplace'){
        x0=(sqrt(mu*mu+4*sd*sd*r)+mu)/2
        x0^r*exp(-(x0-mu)*(x0-mu)/2/sd/sd)*
            sqrt(2*pi/(r/x0/x0+1/sd/sd))*
            (pnorm(upper,x0,1/sqrt(r/x0+1/sd/sd))- pnorm(lower,x0,1/sqrt(r/x0+1/sd/sd))) 
    }else if (approximation=='int2') {
        mTruncNorm.int2(r,mu,sd,lower,upper,...)
    }else if (approximation=='numerical') {
        obj=function(x,r,mu,sd)x^r*exp(-(x-mu)*(x-mu)/2/sd/sd) ## not working for negative values
        n=max(c(length(r), length(mu), length(sd), length(lower), length(upper)))
         mu=rep(as.double(mu), length=n); sd=rep(as.double(sd),length=n); 
        lower=rep(as.double(lower), length=n);     upper=rep(as.double(upper),length=n)
       ans=numeric(n)
        for(i in 1:n) ans[i]=integrate(obj, lower[i], upper[1],r=r[i],mu=mu[i],sd=sd[i])$value  
        ans
    }
}

#dyn.load("intTruncNorm.so")
mTruncNorm.int2=function(r=as.integer(1), mu=0.0, sd=1.0, lower=-Inf, upper=Inf, takeLog=TRUE, ndiv=8)
{   ## return int_lower^upper x^(r-1+1:len) exp(-(x-mu)^2/2/sd^2) dx, ### for integer r only! (using C code)
    n=(max(c(length(r), length(mu), length(sd), length(lower), length(upper))))
    r=rep(r,length=n); 
    mu=rep(as.double(mu), length=n); sd=rep(as.double(sd),length=n); 
    lower=rep(as.double(lower), length=n);     upper=rep(as.double(upper),length=n)
    if(any(lower>upper)) stop("boundary inconsistent")

    integer.r.idx=(r==round(r))
    n.int=as.integer(sum(integer.r.idx))
    ans=numeric(n)
    if(n.int>0)    ans[integer.r.idx]=.C('intTruncNormVec', n=n.int, 
                                               r=as.integer(r[integer.r.idx]), 
                                               mu=mu[integer.r.idx],  sd=sd[integer.r.idx], 
                                               low=lower[integer.r.idx], upp=upper[integer.r.idx], 
                                               ans=ans[integer.r.idx], NAOK=TRUE,PACKAGE='pi0')$ans
    if(n-n.int>0) ans[!integer.r.idx]=.C('fracTruncNormVec', n=n-n.int, 
                                               r=as.double(r[!integer.r.idx]), 
                                               mu=mu[!integer.r.idx],  sd=sd[!integer.r.idx], 
                                               low=lower[!integer.r.idx], upp=upper[!integer.r.idx], 
                                               ans=ans[!integer.r.idx], 
                                               ndiv=as.integer(ndiv[1]), takeLog=as.integer(takeLog[1]), NAOK=TRUE,PACKAGE='pi0')$ans
    return(ans)
}
if(FALSE){
    mTruncNorm(1:3+.2,low=0,approximation='numerical')
    mTruncNorm(1:3+.2,low=0)
}



dt.int2=function(x, df, ncp, log=FALSE, ndiv=8 ) ## pretty fast computation of noncentral t density when df is integer
{   ## when df is integer, this is exact for noncentral t density;
    ## when df is fractional, this is divided difference polynomial interpolation using ndiv points with nearest integer dfs
    if (missing(ncp)) 
        return(dt(x, df, ,log))
	n=max(c(length(x),length(df),length(ncp)))
    x=rep(as.double(x),length=n); df=rep(df,length=n); ncp=rep(as.double(ncp),length=n)

	if(any(finIdx<-is.infinite(df))){
		ans=rep(NA_real_, n)
		ans[finIdx]=dnorm(x[finIdx], ncp[finIdx], , log)
		ans[!finIdx]=Recall(x[!finIdx], df[!finIdx], ncp[!finIdx], log, ndiv)
		return(ans)
	}
    
    integer.df.idx=(df==round(df))
    n.int=as.integer(sum(integer.df.idx))

    df.half=df/2
    tsq.df=x*x+df

    logC=df.half*log(df.half)-0.5*log(pi/2)-lgamma(df.half)-(df.half+.5)*log(tsq.df)-df.half/tsq.df*ncp*ncp
    mus=as.double(x/sqrt(tsq.df)*ncp)

    ints=numeric(n); ints+0

    if(n.int>0)   ints[integer.df.idx]=.C('intTruncNormVec', n.int, 
                                               as.integer(df[integer.df.idx]), 
                                               mus[integer.df.idx],  rep(as.double(1.0), n.int), 
                                               numeric(n.int), rep(Inf,n.int), ans=ints[integer.df.idx], NAOK=TRUE,PACKAGE='pi0')$ans
    if(n-n.int>0) ints[!integer.df.idx]=.C('fracTruncNormVec', n=n-n.int, r=as.double(df)[!integer.df.idx], 
                            mu=mus[!integer.df.idx], sd=rep(as.double(1.0),n-n.int), 
                            low=rep(as.double(0.0), n-n.int), upp=rep(Inf, n-n.int),
                            ans=ints[!integer.df.idx], ndiv=as.integer(ndiv), takeLog=as.integer(1), NAOK=TRUE,PACKAGE='pi0')$ans

    if(log) { logC+log(pmax(0,ints))
    }else exp(logC)*pmax(0,ints)
}

if(FALSE){
 load('work.RData')
 dyn.load("intTruncNorm.so");dyn.unload("intTruncNorm.so");system("R CMD SHLIB intTruncNorm.c");dyn.load("intTruncNorm.so")
 dt.int2(tstat[1],dfs[1],ncps[1])
 dt(tstat[1],dfs[1],ncps[1])
 tmp1=dt.int2(tstat,dfs,ncps)
 tmp0=dt(tstat,dfs,ncps)
 relerr=(tmp1-tmp0)/tmp0
# plot(tmp1,tmp0); abline(0,1,lwd=3,col=3)
plot(dfs, relerr*100); abline(v=round(dfs), col=3, h=0)

system.time(replicate(50, dt.int2(tstat,dfs,ncps)))
system.time(replicate(50, dt(tstat,dfs,ncps)))
}
