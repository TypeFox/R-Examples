grad.hessinv <-
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
    A<-matrix(0, nrow = N+1, ncol = N+1)
    A[1,1]<--n/theta^2+sum(temp^2/(1+(theta-1)*temp)^2)+sum(temp1^2/(1+(theta-1)*temp1)^2)
    A[2:(N+1),1]<--(IN-EN%*%t(p))%*%DQ
    A[1,2:(N+1)]<--t(DQ)
    Omega<-PI%*%solve(IN-(theta-1)^2*PI%*%DD%*%PI)%*%PI
    A22.1<--(Omega+(theta-1)*Omega%*%EN%*%t(DQ)%*%Omega/as.numeric(1-(theta-1)*t(DQ)%*%Omega%*%EN))
    A[1,1]<-1/(A[1,1]-A[1,2:(N+1)]%*%A22.1%*%A[2:(N+1),1])
    A[2:(N+1),2:(N+1)]<-A22.1+A22.1%*%A[2:(N+1),1]%*%A[1,2:(N+1)]%*%A22.1*A[1,1]
    A[1,2:(N+1)]<--A[1,1]*A[1,2:(N+1)]%*%A22.1
    A[2:(N+1),1]<--A22.1%*%A[2:(N+1),1]*A[1,1]
    list(H=c(h0,H), Ainv=A)
}

