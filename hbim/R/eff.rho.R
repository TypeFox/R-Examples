`eff.rho` <-
function(mu,sigma=.553,
    rho=   c(0,      .25,   .5,     .75,   1),
    COLORS=c("black","blue","green","red","black"),...){
    nmu<-length(mu)
    nr<-length(rho)
    out1<-rep(NA,nmu)
    for (i in 1:nmu){
        out1[i]<-1-hbrr(mu[i],sigma^2,...)
    }

    out2<-out3<-matrix(rep(NA,nmu*nr),nmu,nr)
    for (i in 1:nmu){
        for (j in 1:nr){
            V<-make.v(2,rho[j],sigma^2)
            out2[i,j]<-1-hbrr(c(mu[i],mu[i]),V,...)
        }
    }
    for (i in 1:nmu){
        for (j in 1:nr){
            V<-make.v(3,rho[j],sigma^2)
            out3[i,j]<-1-hbrr(c(mu[i],mu[i],mu[i]),V,...)
        }
    }
    out<-list(mu=mu,
         out1=out1,out2=out2,out3=out3,col1=c("black"),col2=COLORS,col3=COLORS,cparms=rho)
    return(out)
}

