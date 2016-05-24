coeff.kzft<-function(m, k){
      poly<-polynomial(rep(1,m))
      polyk<-poly^k
      coef<-as.vector(polyk)
      coef<-coef/m^k
      return(coef)
}

transfun.kzft<-function(m, k, lamda=seq(-0.5,0.5,by=0.01), omega=0){

      lamda<-lamda*2*pi
      omega<-omega*2*pi

      N<-length(lamda)

      tf<-array(0,dim=c(N, m))

      for ( j in (1:m) ){
         tf[,j]<-exp(1i*(lamda-omega)*j)
      }

      tf1<-rep(0,N)
      for ( i in (1:N) ){
         tf1[i]<-sum(tf[i,])
      }

      tf2<-(1/m)*tf1
      tf2<-abs(tf2)^k
      return(tf2)
}

kzft<-function(x,m,k,p=1,n=1){

    data<-as.vector(x)

    if (((m-1)*k+1)>length(data))
        stop("invalid 'm' & 'k':(m-1)k+1 should be less equal length of data")
 
    coef<-coeff.kzft(m,k)

    N<-length(data)
    M<-(m-1)*k+1
    L<-round(M*p)
    T<-floor((N-M)/L)+1
    kzft<-array(NA,dim=c(T,n*m))

    omega<-2*pi*seq(0,1,length=n*m+1)[-1]

    s<-0:(M-1)

    coefft<-coef*exp(-1i*s%o%omega)

    for ( t in (1:T) ){
        kzft[t,]<-data[((t-1)*L+1):((t-1)*L+M)]%*%coefft
    }
    kzftf<-colMeans(kzft)
    lst<-list(tfmatrix=kzft, fft=kzftf)
    return(lst)
}



