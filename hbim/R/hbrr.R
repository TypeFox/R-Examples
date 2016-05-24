`hbrr` <-
function(mu,v,a=rep(1,length(mu)),simulate=FALSE,nsim=10^4,...){
    dim.mu<-length(mu)
    if (dim.mu>3 & !simulate) stop("program only allows up to 3 components for integration")
    if (dim.mu>1){
        if (dim(v)[1]!=dim.mu | dim(v)[2]!=dim.mu) stop("if length(mu)=p, then v should be a p X p matrix")
    }
    if (simulate){
        out<-hbrr.simulate(MU=mu,V=v,A=a,NSIM=nsim)
    }
    else if (dim.mu==1){
        out<-hbrr.integrate1(MU=mu,V=v,A=a,...)
    }
    else if (dim.mu==2){
        if (v[1,2]==v[1,1]) out<-hbrr.integrate2.rhoeq1(MU=mu,V=v,A=a,...)
        else out<-hbrr.integrate2(MU=mu,V=v,A=a,...)
    }
    else if (dim.mu==3){
        if (v[1,2]==v[1,1]) out<-hbrr.integrate3.rhoeq1(MU=mu,V=v,A=a,...)
        else out<-hbrr.integrate3(MU=mu,V=v,A=a,...)
    }
    return(out)
}

