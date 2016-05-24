kzp<-function(x,m,k,p=1,n=1){
    data<-as.vector(x)
    if (((m-1)*k+1)>length(data))
        stop("invalid 'm' & 'k':(m-1)k+1 should be less equal length of data")

    coef<-coeff.kzft(m,k)

    N<-length(data);
    M<-(m-1)*k+1
    L<-round(M*p)    
    T<-floor((N-M)/L)+1

    kzft<-array(NA,dim=c(T,n*m))
    pg<-array(NA,dim=c(T,n*m))


    omega<-2*pi*seq(0,1,length=n*m+1)[-1]

    s<-0:(M-1);

    coefft<-coef*exp(-1i*s%o%omega)

    for ( t in (1:T) ){
        kzft[t,]<-data[((t-1)*L+1):((t-1)*L+M)]%*%coefft
        pg[t,]<-((abs(kzft[t,]))^2)*M
    }

    kzp<-colMeans(pg)

    kzp<-kzp[1:round(n*m/2)]

    return(kzp)
}

nonlinearity.kzp<-function(pg,K=length(pg)){
    N<-length(pg)

    S<-rep(0,N)
    S[1]<-0
    S[N]<-0
    for ( t in 2:(N-1) ) {
        S[t]<-abs(pg[t+1]-2*pg[t]+pg[t-1])
    }

    total<-sum(S)

    sq<-array(0, dim=c(N, K))

    for ( i in (1:N) ) for ( j in (2:K) ) {
        sq[i,j]<-sum(S[(max(1,(i-j+1))):(min(N,(i+j-1)))])
    }

    lst<-list(total=total, sqmatrix=sq)
    return(lst)
}


variation.kzp<-function(pg,K=length(pg)){
    N<-length(pg)
    S=c(diff(pg)^2,0)

    total<-sum(S)

    sq<-array(0, dim=c(N, K))
    for ( i in (1:N) ) for ( j in (2:K) ) {
        sq[i,j]<-sum(S[(max(1,(i-j+1))):(min(N,(i+j-2)))])
    }

    lst<-list(total=total, sqmatrix=sq)
    return(lst)
}

smooth.kzp<-function(pg,c,K=length(pg),method = "DZ") {
    N<-length(pg)
    spg<-rep(0,N)
    m<-rep(0,N)

    if      (method == "DZ") sq<-variation.kzp(pg,K)
    else if (method == "NZ") sq<-nonlinearity.kzp(pg,K)

    cc<-c*sq$total

    for ( i in (1:N) ) {
        m[i]<-sum(sq$sqmatrix[i,1:K]<=cc)
        spg[i]<-mean(pg[(max(1,(i-m[i]+1))):(min(N,(i+m[i]-1)))])
    }

    lst<-list(periodogram=spg,number=m)
    return(lst)
}

