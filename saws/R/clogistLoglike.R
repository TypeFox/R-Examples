`clogistLoglike` <-
function(n,m,x,beta){
    M<-sum(m)
    N<-sum(n)
    if (M==0)  return(0)
    else if (M==N) return(0)
    x<-as.matrix(x)
    eta<- x %*% beta
    U<-exp(eta)
    if (M==1) return(sum(eta*m) - log(sum(U*n)) )
    if (M>N/2){
        ## for efficiency, keep loop part of calculation to minimum
        ## by switching m and n-m, beta and -beta
  m<-n-m
  M<-N-M
  U<-1/U
  eta<- -eta
    }
    if (M==1) return(sum(eta*m) - log(sum(U*n)) )
    B<-rep(1,N-M+1)
    u<-rep(NA,N)
    count<-1
    for (a in 1:length(n)){
        u[count:(count+n[a]-1)]<-U[a]
        count<-count+n[a]
    }
    ## The last 2 lines of this function, may be written more 
    ## clearly (i.e., more like in Gail, et al) BUT LESS EFFICIENTLY as:
    #B<-matrix(0,M+1,N+1)
    #B[1,]<-1
    #for (i in 1:M){
    #for (j in i:(N-M+i)){
    #B[i+1,j+1]<- B[i+1,j]+u[j]*B[i,j]
    #}
    #} 
    #sum(eta*m) - log(B[M+1,N+1])
    for (i in 1:(M-1)) B<- cumsum(B*u[i:(N-M+i)])
    sum(eta*m) - log(sum(B*u[M:N]))
}

