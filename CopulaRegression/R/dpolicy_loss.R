dpolicy_loss <-
function(l,mu,delta,lambda,theta,family,y.max=300,zt=TRUE){
    n<-length(l)
    if (length(lambda)==1) lambda<-rep(lambda,n)
    if (length(mu)==1) mu<-rep(mu,n)
    out<-vector(length=n)
   
    # define copula data for poisson variable
    
    for (i in 1:n){
    if (zt==TRUE){
         y<-1:y.max
	   v<-pztp(y,lambda[i])
        vv<-pztp(y-1,lambda[i])
    }
    if (zt==FALSE){
         y<-1:y.max
        v<-ppois(y,lambda[i])
        vv<-ppois(y-1,lambda[i])
    }
        u<-pgam(l[i]/y,mu[i],delta)
        par_der<-D_u(u,v,theta,family)-D_u(u,vv,theta,family)
        dummy<-par_der*dgam(l[i]/y,mu[i],delta)/y
        out[i]<-sum(dummy)
    }
    out[l<=0]=0
    return(out)
}




