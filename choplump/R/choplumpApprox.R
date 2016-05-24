`choplumpApprox` <-
function(W,Z,use.ranks=TRUE){
    ord<- order(W,Z)
    w<-W[ord]
    z<-Z[ord]
    M<-length(w[w!=0])
    m0<-length(z[w!=0 & z==0])
    m1<-M-m0
    n1<-length(z[z==1])
    N<- length(w)
    K<-N-M
    n0<-N-n1
    k0<-n0-m0
    k1<-n1-m1
    WM<-w[(N+1-M):N]
    if (use.ranks) SM<-rank(WM)
    else SM<-WM
    ZM<-z[(N+1-M):N]

    d<-data.frame(W=w,Z=z)
    dchop<-chopGeneral(d)
    if (use.ranks) T0<-TDiM(rank(dchop$W),dchop$Z)
    else T0<- TDiM(dchop$W,dchop$Z)

    sapply.func<-function(h){
        Qh.times.dhyper(h,n1,n0,M,SM,T0,use.ranks)
    }
    Nperm<-choose(N,n1)
    Ti<-sapply(max(0,n1-M):min(n1,K),sapply.func)
    p.lower<- 1-sum(Ti)
    p.upper<- sum(Ti)
    out<-c(p.lower=p.lower,p.upper=p.upper,p.2sided=min(1,2*min(p.lower,p.upper)))
    out
}

