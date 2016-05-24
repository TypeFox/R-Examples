dFsnc.mix=function(F,df1,df2, delta0, gamma2, log=FALSE, approximation=#c('int2','saddlepoint','laplace','none'),...)
                                                                       'none',...)
{
    approximation=match.arg(approximation)
    if(approximation=='none'){
        scale.fact=(1+gamma2)
        ncp=delta0/scale.fact
        return( if(log){
                    df(F/scale.fact, df1, df2, ncp, log=TRUE)-log(scale.fact)
                }else{
                    df(F/scale.fact, df1, df2, ncp)/(scale.fact)
                }
        )
    }else {
        stop("other approximations not implemented yet.")
    }
#    else if (approximation=='int2'){
#        scale.fact=sqrt(1+sd.ncp*sd.ncp)
#        ncp=mu.ncp/scale.fact
#        return( if(log){
#                    dt.int2(t/scale.fact, df, ncp, log=TRUE,...)-log(scale.fact)
#                }else{
#                    dt.int2(t/scale.fact, df, ncp,...)/(scale.fact)
#                }
#        )
#    }else if (approximation=='saddlepoint'){
#        scale.fact=sqrt(1+sd.ncp*sd.ncp)
#        ncp=mu.ncp/scale.fact
#        return( if(log){
#                    dt.sad(t/scale.fact, df, ncp, log=TRUE,...)-log(scale.fact)
#                }else{
#                    dt.sad(t/scale.fact, df, ncp,...)/(scale.fact)
#                }
#        )
#    }else if (approximation=='laplace'){
#        denom=(1+sd.ncp*sd.ncp)*df+t*t
#        u0=mu.ncp*t*sqrt(df+t*t)/denom
#        g0sq=(df+t*t)*(1+sd.ncp*sd.ncp)/denom
#        g0=sqrt(g0sq)
#        x0=(sqrt(4*g0sq*df+u0*u0)+u0)/2
#
#        exponent=-mu.ncp*mu.ncp*df/2/denom + df*log(x0) -(x0-u0)*(x0-u0)/2/g0sq
#
#        norm.prob=pnorm(0, x0, x0*g0/sqrt(g0sq*df+x0*x0), lower.tail=FALSE, log.p=TRUE)
#        norm.prob0=pnorm(0,g0*sqrt(df),g0/sqrt(2),lower.tail=FALSE,log.p=TRUE)
#
#
#        ans=(exponent
#            +log(x0)+df/2-df/2*log(df+t*t)
#            -lbeta(df/2,.5)
#            -.5*log(.5)-.5*log((1+sd.ncp*sd.ncp)*df+t*t)-.5*log(g0sq*df+x0*x0)
#            +norm.prob - norm.prob0
#        )
#        return (    if(log) ans else exp(ans)   )
#    }
}
