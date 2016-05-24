density_joint <-
function(x,y,mu,delta,lambda,theta,family,zt=TRUE){
    if (length(mu)!=length(lambda)) stop("mu and lambda must be of the same length")
    u<-pgam(x,mu,delta)
    u[u>=1]=1
    u[u<0]=0
    u[is.na(u)]=0
    u[u==-Inf]=0
    u[u==Inf]=1
    if (zt==TRUE){
        v<-pztp(y,lambda)
        vv<-pztp(y-1,lambda)
    }
    if (zt==FALSE){
        v<-ppois(y,lambda)
        vv<-ppois(y-1,lambda)
    }
    v[v>=1]=1
    v[v<0]=0
    v[is.na(v)]=0
    v[v==-Inf]=0
    v[v==Inf]=1
    vv[vv>=1]=1
    vv[vv<0]=0
    vv[is.na(vv)]=0
    vv[vv==-Inf]=0
    vv[vv==Inf]=1
    marginal.x<-dgam(x,mu,delta)
    # compute partial derivative
    par_der<-D_u(u,v,theta,family)
    par_der1=D_u(u,vv,theta,family)
    dummy<-par_der - par_der1
    if (zt==TRUE){
        dummy[y==1]=par_der[y==1]
    }
    if (zt==FALSE){
        dummy[y==0]=par_der[y==0]
    }
    # multiply with density of X
    out<-marginal.x*dummy
    return(out)
}
