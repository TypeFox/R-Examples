`eff.mu` <-
function(mu,factor=c(1/10,1/3,1/2,1),
    COLORS=c("red","green","blue","black"),
    sigma=.553,rho=0,...){
    nmu<-length(mu)
    out1<-rep(NA,nmu)
    for (i in 1:nmu){
        out1[i]<-1-hbrr(mu[i],sigma^2,...)
    }
    nf<-length(factor)
    
    out2<-out3<-mu2<-matrix(rep(NA,nmu*nf),nmu,nf)
    V<-make.v(2,rho,sigma^2)
    
    for (j in 1:nf){
        mu2[,j]<- mu+log10(factor[j])
        for (i in 1:nmu){
            out2[i,j]<-1-hbrr(c(mu[i],mu2[i,j]),V,...)
        }
    }
    V<-make.v(3,rho,sigma^2)
    for (j in 1:nf){
        for (i in 1:nmu){
            out3[i,j]<-1-hbrr(c(mu[i],mu2[i,j],mu2[i,j]),V,...)
        }
    }
    out<-list(mu=mu,mu2=mu2,
         out1=out1,out2=out2,out3=out3,col1=c("black"),col2=COLORS,col3=COLORS,cparms=factor,sigma=sigma,rho=rho)
    return(out)
}

