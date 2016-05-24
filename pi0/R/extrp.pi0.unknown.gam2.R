extrp.pi0.gam2=function(n1,n2,y,gam2.interval=c(1e-3,6))
{   ## y=(1-pi0)*B+pi0  <==> y-B=(1-B)*pi0
    f.obj=function(parm){gam2=parm
        B=sqrt((n1+n2)/(n1+n2+n1*n2*gam2))
        pi0.ans=min(1,max(0,crossprod(y-B,1-B)/crossprod(1-B)))
        sse.ans=crossprod(y-(1-pi0.ans)*B-pi0.ans)
        sse.ans
    }
    optimize.fit=optimize(f.obj,gam2.interval)
    sse.ans=optimize.fit$objective
    gam2=optimize.fit$minimum
    B=sqrt((n1+n2)/(n1+n2+n1*n2*gam2))
    pi0.ans=min(1,max(0,crossprod(y-B,1-B)/crossprod(1-B)))

    list(
        pi0=pi0.ans,
        sse=sse.ans, 
        gamma2=gam2,
        par=c(pi0=pi0.ans,gamma2=gam2,slope=1-pi0.ans,rate=0.5),
        fitted.obj=optimize.fit
    )
}

extrp.pi0.slope.gam2=function(n1,n2,y,gam2.interval=c(1e-3,6),eps=1e-5)
{
    #library("limSolve")
    my=mean(y)
    Gmat=matrix(c(1,0,
                 -1,0,
                  0,1),3,2,byrow=TRUE)
    hvec=matrix(c(0,-1,eps))
    
    f.obj=function(parm){gam2=parm
        xvec=sqrt((n1+n2)/(n1+n2+n1*n2*gam2))
        mx=mean(xvec);
        slope.ls=crossprod(xvec-mx,y-my)/crossprod(xvec-mx)
        pi0.ls=my-slope.ls*mx
        if(slope.ls>0 && pi0.ls >=0 && pi0.ls<=1)
            return(crossprod(y-slope.ls*xvec-pi0.ls)
#                   list(pi0=pi0.ls,
#                        slope=slope.ls,
#                        sse=crossprod(y-slope.ls*xvec-pi0.ls),
#                        par=c(pi0=pi0.ls,gamma2=gam2,slope=slope.ls,rate=0.5),
#                        fitted.obj=NULL)
            )

        ## QP solution with constraints
    #    cat("LS=",pi0.ls,slope.ls,crossprod(y-slope.ls*xvec-pi0.ls),fill=T)
        X=cbind(1,c(xvec))
        lsei.fit=limSolve::lsei(A=X,B=y,G=Gmat,H=hvec,verbose=FALSE)
        return(lsei.fit$solutionNorm
#               list(pi0=lsei.fit$X[1],
#                    slope=lsei.fit$X[2],
#                    sse=lsei.fit$solutionNorm,
#                    par=c(pi0=lsei.fit$X[1],gamma2=gam2,slope=lsei.fit$X[2],rate=0.5),
#                    fitted.obj=NULL)
        )
    }

    optimize.fit=optimize(f.obj,gam2.interval)
    sse.ans=optimize.fit$objective
    gam2=optimize.fit$minimum

    xvec=sqrt((n1+n2)/(n1+n2+n1*n2*gam2))
    mx=mean(xvec);
    slope.ls=crossprod(xvec-mx,y-my)/crossprod(xvec-mx)
    pi0.ls=my-slope.ls*mx
    if(slope.ls>0 && pi0.ls >=0 && pi0.ls<=1)
        return(list(pi0=pi0.ls,
                    slope=slope.ls,
                    gamma2=gam2,
                    sse=crossprod(y-slope.ls*xvec-pi0.ls),
                    par=c(pi0=pi0.ls,gamma2=gam2,slope=slope.ls,rate=0.5),
                    fitted.obj=NULL)
        )
    X=cbind(1,c(xvec))
    lsei.fit=limSolve::lsei(A=X,B=y,G=Gmat,H=hvec,verbose=FALSE)
    return(list(pi0=lsei.fit$X[1],
                slope=lsei.fit$X[2],
                gamma2=gam2,
                sse=lsei.fit$solutionNorm,
                par=c(pi0=lsei.fit$X[1],gamma2=gam2,slope=lsei.fit$X[2],rate=0.5),
                fitted.obj=NULL)
    )

}

extrp.pi0.rate.gam2=function(n1,n2,y,gam2.interval=c(1e-3,6),rate.interval=c(.3,2),eps=1e-5)
{   ## y=(1-pi0)*B+pi0  <==> y-B=(1-B)*pi0
    f.obj0=function(parm0,parm1){#       rate=parm0
        gam2=parm1
        B=((n1+n2)/(n1+n2+n1*n2*gam2))^parm0
        pi0.ans=min(1,max(0,crossprod(y-B,1-B)/crossprod(1-B)))
        sse.ans=crossprod(y-(1-pi0.ans)*B-pi0.ans)
        sse.ans
    }
    f.obj=function(parm){#    gam2=parm=parm1
        optimize.fit0=optimize(f.obj0,rate.interval,parm1=parm)
        sse=optimize.fit0$objective
        rate=optimize.fit0$minimum
        attr(sse,'rate')=rate
        sse
    }

    optimize.fit=optimize(f.obj,gam2.interval)
    sse.ans=optimize.fit$objective
    gam2=optimize.fit$minimum
    rate=attr(sse.ans,'rate')


    B=((n1+n2)/(n1+n2+n1*n2*gam2))^rate
    pi0.ans=min(1,max(0,crossprod(y-B,1-B)/crossprod(1-B)))
    return(list(
            pi0=pi0.ans,
            rate=rate,
            gamma2=gam2,
            sse=sse.ans,
            par=c(pi0=pi0.ans,gamma2=gam2,slope=1-pi0.ans,rate=rate),
            fitted.obj=optimize.fit)
    )
}

extrp.pi0.both.gam2=function(n1,n2,y,gam2.interval=c(1e-3,6),rate.interval=c(.3,2),eps=1e-5)
{   ## y=(1-pi0)*B+pi0  <==> y-B=(1-B)*pi0
        #library("limSolve")
    my=mean(y)
    Gmat=matrix(c(1,0,
                 -1,0,
                  0,1),3,2,byrow=TRUE)
    hvec=matrix(c(0,-1,eps))
    
    f.obj0=function(parm0,parm1){#       rate=parm0
        gam2=parm1
        xvec=((n1+n2)/(n1+n2+n1*n2*gam2))^parm0

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

    f.obj=function(parm){#    gam2=parm=parm1
        optimize.fit0=optimize(f.obj0,rate.interval,parm1=parm)
        sse=optimize.fit0$objective
        rate=optimize.fit0$minimum
        attr(sse,'rate')=rate
        sse
    }

    optimize.fit=optimize(f.obj,gam2.interval)
    sse.ans=optimize.fit$objective
    gam2=optimize.fit$minimum
    rate=attr(sse.ans,'rate')

    ## last fit
    xvec=((n1+n2)/(n1+n2+n1*n2*gam2))^rate
    mx=mean(xvec);
    slope.ls=crossprod(xvec-mx,y-my)/crossprod(xvec-mx)
    pi0.ls=my-slope.ls*mx
    if(slope.ls>0 && pi0.ls >=0 && pi0.ls<=1)
        return(list(pi0=pi0.ls,
                    slope=slope.ls,
                    rate=rate,
                    gamma2=gam2,
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
