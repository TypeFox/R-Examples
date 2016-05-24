density_conditional <-
function(y,x,mu,delta,lambda,theta,family,zt=TRUE){
    # define copula data
    x<-rep(x,length(y))
    u=pgam(x,mu,delta)
    if (zt==TRUE){
        v=pztp(y,lambda)
        vv<-pztp(y-1,lambda)
    }
    if (zt==FALSE){
        v=ppois(y,lambda)
        vv<-ppois(y-1,lambda)
    }
    #cat(paste("v\n"))
    #cat(paste(v,"\n"))
    #cat(paste("vv\n"))
    #cat(paste(vv,"\n"))
    # compute partial derivative
    out<-D_u(u,v,theta,family)- D_u(u,vv,theta,family) 
    if (zt==TRUE){
	out[y==1]=(D_u(u,v,theta,family))[y==1]
    }
    if (zt==FALSE){
	out[y==0]=(D_u(u,v,theta,family))[y==0]
    }
    out[out<0]=0
    out[out>1]=1
    return(out)

}
