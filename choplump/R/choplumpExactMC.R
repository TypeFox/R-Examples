`choplumpExactMC` <-
function(W,Z,use.ranks=TRUE,nMC=10^4-1,seed=1234321){
    set.seed(seed)
    ord<- order(W,Z)
    w<-W[ord]
    z<-Z[ord]
    M<-length(w[w!=0])
    m0<-length(z[w!=0 & z==0])
    m1<-M-m0
    n1<-length(z[z==1])
    N<- length(w)
    n0<-N-n1
    k0<-n0-m0
    k1<-n1-m1
    WM<-w[(N+1-M):N]
    if (use.ranks){
        RM<-rank(WM)
        NTIES<-table(RM)
        testfunc<-function(ZM){
            testfunc.wilcox.ties(chop(ZM,n1,n0,M),RM,NTIES,M)
        }
    }
    else {
        testfunc<-function(ZM){
            testfunc.DiM(chop(ZM,n1,n0,M),WM)
        }
   }
    T0<-testfunc(z[(N+1-M):N]) 
    Ti<-rep(NA,nMC)

    for (i in 1:nMC){
        Ti[i]<-testfunc(sample(z)[(N+1-M):N])
    }
    S.lte <- length((1:nMC)[Ti <= T0])
    S.gte <- length((1:nMC)[Ti >= T0])
    p.lower <- (S.lte + 1)/(nMC + 1)
    p.upper <- (S.gte + 1)/(nMC + 1)
    out<-c(p.lower=p.lower,p.upper=p.upper,p.2sided=min(1,2*min(p.lower,p.upper)))
    out
}

