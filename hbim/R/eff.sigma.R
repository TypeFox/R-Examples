`eff.sigma` <-
function(mu,sigmas,COLORS=c("red","green","blue"),rho=0,...){
    nmu<-length(mu)
    nsig<-length(sigmas)
    out1<-out2<-out3<-matrix(rep(NA,nmu*nsig),nmu,nsig)
    for (i in 1:nmu){
        for (j in 1:nsig){
            out1[i,j]<-1-hbrr(mu[i],sigmas[j]^2,...)
        } 
    }

    for (i in 1:nmu){
        for (j in 1:nsig){
            V<-make.v(2,rho,sigmas[j]^2)
            out2[i,j]<-1-hbrr(c(mu[i],mu[i]),V,...)
        }
    }

    for (i in 1:nmu){
        for (j in 1:nsig){
            V<-make.v(3,rho,sigmas[j]^2)
            out3[i,j]<-1-hbrr(c(mu[i],mu[i],mu[i]),V,...)
        } 
    }
    out<-list(mu=mu,
         out1=out1,out2=out2,out3=out3,col1=COLORS,col2=COLORS,col3=COLORS,cparms=sigmas)
    return(out)
}

