extrp.pi0.only=function(n1,n2,y,gam2)
{   ## y=(1-pi0)*B+pi0  <==> y-B=(1-B)*pi0
    B=sqrt((n1+n2)/(n1+n2+n1*n2*gam2))
    pi0.ans=min(1,max(0,crossprod(y-B,1-B)/crossprod(1-B)))
    sse.ans=crossprod(y-(1-pi0.ans)*B-pi0.ans)
    list(
        pi0=pi0.ans,
        sse=sse.ans, 
        par=c(pi0=pi0.ans,gamma2=gam2,slope=1-pi0.ans,rate=0.5),
        fitted.obj=NULL
    )
}

extrp.pi0.slope=function(n1,n2,y,gam2,eps=1e-5)
{
    xvec=sqrt((n1+n2)/(n1+n2+n1*n2*gam2))
    mx=mean(xvec);my=mean(y)
    slope.ls=crossprod(xvec-mx,y-my)/crossprod(xvec-mx)
    pi0.ls=my-slope.ls*mx
    if(slope.ls>0 && pi0.ls >=0 && pi0.ls<=1)
        return(list(pi0=pi0.ls,
                    slope=slope.ls,
                    sse=crossprod(y-slope.ls*xvec-pi0.ls),
                    par=c(pi0=pi0.ls,gamma2=gam2,slope=slope.ls,rate=0.5),
                    fitted.obj=NULL))

    ## QP solution with constraints
    #library("limSolve")
#    cat("LS=",pi0.ls,slope.ls,crossprod(y-slope.ls*xvec-pi0.ls),fill=T)
    X=cbind(1,c(xvec))
    Gmat=matrix(c(1,0,
                 -1,0,
                  0,1),3,2,byrow=TRUE)
    hvec=matrix(c(0,-1,eps))
    lsei.fit=limSolve::lsei(A=X,B=y,G=Gmat,H=hvec,verbose=FALSE)
    return(list(pi0=lsei.fit$X[1],
                slope=lsei.fit$X[2],
                sse=lsei.fit$solutionNorm,
                par=c(pi0=lsei.fit$X[1],gamma2=gam2,slope=lsei.fit$X[2],rate=0.5),
                fitted.obj=NULL)
    )
}

extrp.pi0.rate=function(n1,n2,y,gam2,rate.interval=c(.3,2),eps=1e-5)
{   ## y=(1-pi0)*B+pi0  <==> y-B=(1-B)*pi0
    f.obj=function(parm){#       rate=parm
        B=((n1+n2)/(n1+n2+n1*n2*gam2))^parm
        pi0.ans=min(1,max(0,crossprod(y-B,1-B)/crossprod(1-B)))
        sse.ans=crossprod(y-(1-pi0.ans)*B-pi0.ans)
        sse.ans
    }

    optimize.fit=optimize(f.obj,rate.interval)
    sse=optimize.fit$objective
    rate=optimize.fit$minimum

    B=((n1+n2)/(n1+n2+n1*n2*gam2))^rate
    pi0.ans=min(1,max(0,crossprod(y-B,1-B)/crossprod(1-B)))
    return(list(
            pi0=pi0.ans,
            rate=rate,
            sse=sse,
            par=c(pi0=pi0.ans,gamma2=gam2,slope=1-pi0.ans,rate=rate),
            fitted.obj=optimize.fit)
    )
}

extrp.pi0.both=function(n1,n2,y,gam2,rate.interval=c(.3,2),eps=1e-5)
{   ## y=(1-pi0)*B+pi0  <==> y-B=(1-B)*pi0
        #library("limSolve")
    my=mean(y)
    Gmat=matrix(c(1,0,
                 -1,0,
                  0,1),3,2,byrow=TRUE)
    hvec=matrix(c(0,-1,eps))
    f.obj=function(parm){#       rate=parm
        xvec=((n1+n2)/(n1+n2+n1*n2*gam2))^parm

        mx=mean(xvec);
        slope.ls=crossprod(xvec-mx,y-my)/crossprod(xvec-mx)
        pi0.ls=my-slope.ls*mx
        if(slope.ls>0 && pi0.ls >=0 && pi0.ls<=1)
            return(crossprod(y-slope.ls*xvec-pi0.ls))

        ## QP solution with constraints
        X=cbind(1,c(xvec))
        lsei.fit=limSolve::lsei(A=X,B=y,G=Gmat,H=hvec,verbose=FALSE)
        return(lsei.fit$solutionNorm)
    }

    optimize.fit=optimize(f.obj,rate.interval)
    sse=optimize.fit$objective
    rate=optimize.fit$minimum

    ## last fit
    xvec=((n1+n2)/(n1+n2+n1*n2*gam2))^rate
    mx=mean(xvec);
    slope.ls=crossprod(xvec-mx,y-my)/crossprod(xvec-mx)
    pi0.ls=my-slope.ls*mx
    if(slope.ls>0 && pi0.ls >=0 && pi0.ls<=1)
        return(list(pi0=pi0.ls,
                    slope=slope.ls,
                    rate=rate,
                    sse=crossprod(y-slope.ls*xvec-pi0.ls),
                    par=c(pi0=pi0.ls,gamma2=gam2,slope=slope.ls,rate=rate),
                    fitted.obj=optimize.fit))
    X=cbind(1,c(xvec))
    lsei.fit=limSolve::lsei(A=X,B=y,G=Gmat,H=hvec,verbose=FALSE)
    return(list(pi0=lsei.fit$X[1],
                slope=lsei.fit$X[2],
                rate=rate,
                sse=lsei.fit$solutionNorm,
                par=c(pi0=lsei.fit$X[1],gamma2=gam2,slope=lsei.fit$X[2],rate=rate),
                fitted.obj=optimize.fit)
    )
}
