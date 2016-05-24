
dt.sad=function(x, df, ncp=0, log=FALSE, normalize=c('approximate','derivative','integral','none'),epsilon=1e-4){
    if(missing(df)) stop('Error: the number of degrees of freedom must be supplied!')
    normalize=match.arg(normalize)

    unq.dfmu=unique(cbind(df,ncp))
    N=max(c(length(x),length(df),length(ncp)))
    y1=rep(x,length=N)
    mu=ncp
    n=df

    y2.hat=(mu*y1+sqrt(4*n*(y1*y1+n)+mu*mu*y1*y1))/2/(y1*y1+n)
    t1.hat=-mu+y1*y2.hat
#    t2.hat=-y1*t1.hat/2/n/y2.hat
    w=sqrt(-mu*t1.hat-2*n*log(y2.hat))*ifelse(y1>mu,1,-1)
    u=sqrt(mu*y1*y2.hat/2/n+1)/y2.hat

    ans=dnorm(w)/u
    if (normalize=='none') {
        if(log)return(log(ans)) else return(ans)
    }else if(normalize=='approximate') {
        fact=numeric(nrow(unq.dfmu))

        for(i in 1:nrow(unq.dfmu)){
            fact[i]=dt.sad(unq.dfmu[i,2],unq.dfmu[i,1],unq.dfmu[i,2],log=FALSE,normalize='none')/
                        dt(unq.dfmu[i,2],unq.dfmu[i,1],unq.dfmu[i,2])
        }

        tmp1=outer(unq.dfmu[,1],rep(n,length=N),'==')
        tmp2=outer(unq.dfmu[,2],rep(mu,length=N),'==')
        idx=matrix(1:nrow(unq.dfmu),nrow(unq.dfmu),N)[tmp1&tmp2]
        ans=ans/fact[idx]
        if(log)return(log(ans)) else return(ans)
    }else if (normalize=='integral') { 

        fact=numeric(nrow(unq.dfmu))
        for(i in 1:nrow(unq.dfmu)){
            fact[i]=integrate(dt.sad,-Inf,Inf,df=unq.dfmu[i,1],ncp=unq.dfmu[i,2],normalize='none')$value
        }

        tmp1=outer(unq.dfmu[,1],rep(n,length=N),'==')
        tmp2=outer(unq.dfmu[,2],rep(mu,length=N),'==')
        idx=matrix(1:nrow(unq.dfmu),nrow(unq.dfmu),N)[tmp1&tmp2]
        ans=ans/fact[idx]

        if(log)return(log(ans)) else return(ans)
    }else if (normalize=='derivative'){ ## not working properly yet
        d=1/t1.hat/y2.hat      
##########
#        nu=1/(1-2*t2.hat)
#        c1=2*n*nu*nu
#        c2=8*n*nu*nu*nu
#        c3=t1.hat+y1*y2.hat
#        c4=2*n*t2.hat+y1*y1+4*n*n*y2.hat*y2.hat/c1
#        c5=c3/c4
#        c6=2*n/c1*c2*c5*y2.hat*(2*n*t2.hat+y1*y1)+12*n*n*c5*y2.hat+2*c1*y1
#        d.prime=-d*d*(y2.hat*y2.hat-c3*c5)
#        u.prime=-u*c5*(c6/2/c1/c3+2/y2.hat)
######### the above are from the paper but do not work correctly

######### Below are from the direct differentiation in maxima, which works properly
        y1.2=y1*y1
        y2.hat2=y2.hat*y2.hat
        mu2=mu*mu
        A2=(4*n+mu2)*y1.2+4*n*n
        A=sqrt(A2)

        dy2.hat=-(mu*(y1.2-n)*A+(A2-mu2*n)*y1)/2/A/(y1.2+n)/(y1.2+n)
        u.prime=(mu*y2.hat2-(mu*y1*y2.hat+4*n)*dy2.hat)/2/y2.hat2/sqrt(2*n*(mu*y1*y2.hat+2*n))
        dt1.hat=y1*dy2.hat+y2.hat
        d.prime=-dy2.hat/t1.hat/y2.hat/y2.hat-dt1.hat/t1.hat/t1.hat/y2.hat
########## 

        ans=dnorm(w)*(1/u-1/d/w/w/w-d.prime/u+d*u.prime/u/u)

        alpha=mu
        idx=abs(y1-alpha)<epsilon
        if(any(idx)){
            n2=n*n
            n3=n2*n
            n4=n2*n2
            n5=n2*n3
            o1=6*n*n-n
            o2=144*n*n-33*n
            o3=4*n*n-n
            mu4=mu2*mu2
            mu6=mu2*mu4
            lim.F.prime=(2*mu6*n*o1+12*mu4*o1*n2+mu2*n3*o2+24*n4*o3)/(12*sqrt(pi*n5*(mu2+2*n)^7))

            ans[idx]=lim.F.prime[idx]
        }
        if(log)return(log(ans)) else return(ans)
    }
}

if(FALSE){
#curve(dt(x,3,2),-2,15)
#curve(dt.sad(x,3,2,normalize='none'),-2,15,add=T,lty=2)
#curve(dt.sad(x,3,2,normalize='approximate'),-2,15,add=T,col=2)
#curve(dt.sad(x,3,2,normalize='integral'),-2,15,add=T,col=3)
#curve(dt.sad(x,3,2,normalize='derivative'),-2,15,add=T,col=4)
curve(dt.sad(x,1,2,normalize='none'),-2,15,lty=2)
curve(dt(x,1,2),-2,15,add=T)
curve(dt.sad(x,1,2,normalize='approximate'),-2,15,add=T,col=2)
curve(dt.sad(x,1,2,normalize='integral'),-2,15,add=T,col=3)
curve(dt.sad(x,1,2,normalize='derivative'),-2,15,add=T,col=4)
#gc();system.time(dt.sad(seq(-20,20,length=5e5),3,2,normalize='none'))
#gc();system.time(dt.sad(seq(-20,20,length=5e5),3,2,normalize='approximate'))
#gc();system.time(dt.sad(seq(-20,20,length=5e5),3,2,normalize='integral'))
#gc();system.time(dt.sad(seq(-20,20,length=5e5),3,2,normalize='derivative'))
#gc();system.time(dt(seq(-20,20,length=5e5),3,2))
}


pt.sad=function(q, df, ncp=0, log=FALSE,epsilon=1e-4){
    if(missing(df)) stop('Error: the number of degrees of freedom must be supplied!')

    N=max(c(length(q),length(df),length(ncp)))
    y1=rep(q,length=N)
    mu=ncp
    n=n

    y2.hat=(mu*y1+sqrt(4*n*(y1*y1+n)+mu*mu*y1*y1))/2/(y1*y1+n)
    t1.hat=-mu+y1*y2.hat
    t2.hat=-y1*t1.hat/2/n/y2.hat
    w=sqrt(-mu*t1.hat-2*n*log(y2.hat))*ifelse(y1>mu,1,-1)
    u=sqrt(mu*y1*y2.hat/2/n+1)/y2.hat
    d=1/t1.hat/y2.hat
    alpha=mu

    ans=pnorm(w)+dnorm(w)*(1/w-d/u)

    idx=abs(y1-alpha)<epsilon
    if(any(idx)){
        lim.F=1/2-1/6/sqrt(n*pi*(mu*mu+2*n)^3)*mu*(2*mu*mu+3*n)
        n2=n*n
        n3=n2*n
        n4=n2*n2
        n5=n2*n3
        o1=6*n*n-n
        o2=144*n*n-33*n
        o3=4*n*n-n
        mu2=mu*mu
        mu4=mu2*mu2
        mu6=mu2*mu4
        lim.F.prime=(2*mu6*n*o1+12*mu4*o1*n2+mu2*n3*o2+24*n4*o3)/(12*sqrt(pi*n5*(mu2+2*n)^7))

        ans[idx]=(lim.F+(y1-alpha)*lim.F.prime)[idx]
    }
    if(log)return(log(ans)) else return(ans)
}


if(FALSE){
#curve(pt(x,3,2),-2,15,lwd=2)
#curve(pt.sad(x,3,2),-2,15,add=T,col=2)
#
#curve(pt(x,3,3),2.9995,3.0005)
#curve(pt.sad(x,3,3),2.999,3.001,1000,add=T,col=4)
#
gc();system.time(pt.sad(seq(-20,20,length=5e5),3,2))
gc();system.time(pt(seq(-20,20,length=5e5),3,2))
}


