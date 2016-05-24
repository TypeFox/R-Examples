H.Binv <-
function(theta, p, r){
    n<-length(r)
    N<-length(p)
    PI<-diag(as.vector(p))
    temp<-NULL
    temp1<-NULL
    for(j in 1:n){
        temp[j]<-sum(p[1:r[j]])
        temp1[j]<-ifelse(r[j]==1,0,sum(p[1:(r[j]-1)]))
    }
    En<-rep(1,n)
    EN<-rep(1,N); IN<-diag(N)
    D<-matrix(0, nrow = N, ncol = n)
    Q<-1/(1+(theta-1)*temp)
    D1<-matrix(0, nrow = N, ncol = n)
    Q1<-1/(1+(theta-1)*temp1)
    for(k in 1:N) {
        D[k,]<-(r>=k)/(1+(theta-1)*temp)
        D1[k,]<-(r>k)/(1+(theta-1)*temp1)
    }
    DE<-(D+D1)%*%En; DQ<-D%*%Q+D1%*%Q1; DD<-D%*%t(D)+D1%*%t(D1)
    h0<-t(p)%*%DE
    H<-1/p-N-(theta-1)*(DE-h0*EN)
    h0<-n/theta-h0
    Omega<-PI%*%solve(IN-(theta-1)^2*PI%*%DD%*%PI)%*%PI
    Binv<--(Omega+(theta-1)*Omega%*%EN%*%t(DQ)%*%Omega/as.numeric(1-(theta-1)*t(DQ)%*%Omega%*%EN))
    ell<-n*log(theta)+sum(log(p))-sum(log(1+(theta-1)*temp))-sum(log(1+(theta-1)*temp1))
    Fisher<--n/theta^2+t(p)%*%DD%*%p 
    Fisher<-Fisher-t(DQ)%*%Binv%*%(IN-EN%*%t(p))%*%DE
    list(H=H, Binv=Binv, ell=ell)#, Fisher=-Fisher/N)
}
