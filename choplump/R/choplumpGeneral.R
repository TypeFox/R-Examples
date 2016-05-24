`choplumpGeneral` <-
function(W,Z,testfunc=testfunc.wilcox.ties.general){
    w<-W
    z<-Z
    M<-length(w[w!=0])
    m0<-length(z[w!=0 & z==0])
    m1<-M-m0
    n1<-length(z[z==1])
    N<- length(w)
    n0<-N-n1
    k0<-n0-m0
    k1<-n1-m1
    d<-data.frame(W=w,Z=z)
    T0<-testfunc(chopGeneral(d))
    cm<-chooseMatrix(N,n1)
    Nperm<-dim(cm)[1]
    Ti<-rep(NA,Nperm)
    for (i in 1:Nperm){
        d$Z<-cm[i,]
        Ti[i]<-testfunc(chopGeneral(d))
       #Ti[i]<-testfunc()
    }
    p.lower<- length(Ti[Ti<=T0])
    p.upper<- length(Ti[Ti>=T0])
    out<-c(p.lower=p.lower/Nperm,p.upper=p.upper/Nperm,p.2sided=min(1,2*min(p.lower,p.upper)/Nperm))
    out
}

