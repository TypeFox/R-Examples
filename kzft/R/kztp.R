kztp<-function(x,m,k,p=1,n=1,rp1=0,rp2=0.5,cp1=0,cp2=0.5){
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

    rm1<-max(round(m*n*rp1),1)
    rm2<-round(m*n*rp2)
    cm1<-max(round(m*n*cp1),1)
    cm2<-round(m*n*cp2)
    delta.rm<-rm2-rm1+1
    delta.cm<-cm2-cm1+1

    kztp<-array(NA,dim=c(delta.rm,delta.cm,T))
    for ( i in (1:delta.rm) ) for ( j in (1:delta.cm) ){
         kztp[i,j,]<-kzft[,i+rm1-1]*kzft[,j+cm1-1]*Conj(kzft[,i+j+rm1+cm1-2])*M^2
    }               
        
    kztp<-rowMeans(kztp,dim=2)

    lst<-list(bispectrum=kztp, modulus=Mod(kztp), argument=Arg(kztp))
    return(lst)
}

variation.kztp<-function(pg, K=dim(pg)[1]){
    N<-dim(pg)[1]
    sq<-array(0,dim=c(N,N,K))

    sqvar2<-function(pg){
        sdif1<-sum(diff(pg)^2)
        sdif2<-sum(diff(t(pg))^2)

        sumsq<-sdif1+sdif2
        return(sumsq)
    }

    total<-sqvar2(pg)

    for ( i in (1:N) ) for ( j in (1:N) ) for ( k in (2:K) ) {
        sq[i,j,k]<-sqvar2(pg[(max(1,(i-k+1)):min(N,(i+k-1))),(max(1,(j-k+1)):min(N,(j+k-1)))])
    }

    lst<-list(total=total, sqmatrix=sq)
    return(lst)
}

smooth.kztp<-function(pg,c,K=dim(pg)[1]) {
    N<-dim(pg)[1]
    spg<-array(0,dim=c(N,N))
    m<-array(0,dim=c(N,N))

    sq<-variation.kztp(pg,K)

    cc<-c*sq$total


    for ( i in (1:N) ) for ( j in (1:N) ) {
        m[i,j]<-sum(sq$sqmatrix[i,j,1:K]<=cc)
        spg[i,j]<-mean(pg[(max(1,(i-m[i,j]+1))):(min(N,(i+m[i,j]-1))),(max(1,(j-m[i,j]+1))):(min(N,(j+m[i,j]-1)))])
    }

    lst<-list(bispectrum=spg,number=m)
    return(lst)
}


