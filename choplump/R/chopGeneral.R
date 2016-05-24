`chopGeneral` <-
function(d, W="W", Z="Z"){
    ord<- order(d[,W],d[,Z])
    w<-d[ord,W]
    z<-d[ord,Z]
    M<-length(w[w!=0])
    m0<-length(z[w!=0 & z==0])
    m1<-M-m0
    n1<-length(z[z==1])
    N<- length(w)
    n0<-N-n1
    k0<-n0-m0
    k1<-n1-m1

    if (m0/n0 >= m1/n1){
        a<-0
        b<-k1 - floor((n1*k0)/n0 )
    }
    else {
        a<-k0 - floor( (n0*k1)/n1)
        b<-0
    }
    wout<-w[(N+1-(M+a+b)):N]
    zout<-c(rep(0,a),rep(1,b),z[(N+1-M):N])
    dout<-data.frame(W=wout,Z=zout)
    dout
}

