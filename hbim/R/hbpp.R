`hbpp` <-
function(mu,v,a=rep(1,length(mu)),rp=.1,simulate=FALSE,nsim=10^5,...){
    dim.mu<-length(mu)
    if (dim.mu>3 & !simulate) stop("program only allows up to 3 components for integration")
    if (dim.mu>1){
        if (dim(v)[1]!=dim.mu | dim(v)[2]!=dim.mu) stop("if length(mu)=p, then v should be a p X p matrix")
    }
    if (simulate){
        out<-hbpp.simulate(MU=mu,V=v,A=a,RP=rp,NSIM=nsim)
    }
    else if (dim.mu==1){
        out<-hbpp.integrate1(MU=mu,V=v,A=a,RP=rp,...)
    }
    else if (dim.mu==2){
        out<-hbpp.integrate2(MU=mu,V=v,A=a,RP=rp,...)
    }
    else if (dim.mu==3){
        out<-hbpp.integrate3(MU=mu,V=v,A=a,RP=rp,...)
    }
    ## output in percent 
    out<-out*100
    return(out)
}

