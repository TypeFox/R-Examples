`chop` <-
function(ZM,in.n1,in.n0,in.M){
    m0<-length(ZM[ZM==0])
    m1<-in.M-m0
    k0<-in.n0-m0
    k1<-in.n1-m1
    if (m0/in.n0 >= m1/in.n1){
        a<-0
        b<-k1 - floor((in.n1*k0)/in.n0 )
    }
    else {
        a<-k0 - floor( (in.n0*k1)/in.n1)
        b<-0
    }
    out<-list(a=a,b=b,ZM=ZM)
    out
}

