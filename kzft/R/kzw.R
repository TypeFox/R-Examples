
kzw<-function(x,f1=1/length(x),f2=0.5,delta.f=1/length(x),t1=1,t2=length(x),delta.t=1,n,k=3,method="zero"){

updown<-function(x){
    n<-length(x)
    y<-rep(NA,n)

    for ( i in (1:n)) {
     y[i] = x[n-i+1]
    }
    return(y)
}

    data<-as.vector(x)

    left.extra<-floor((round(n/f1)-1)*k/2)
    right.extra<-ceiling((round(n/f1)-1)*k/2)

    left.n.extra<-floor(left.extra/length(data)/2)+1
    right.n.extra<-floor(right.extra/length(data)/2)+1    

    data.left<-rep(c(data,updown(data)),left.n.extra)
    data.right<-rep(c(updown(data),data),right.n.extra)
    
    data_left <-updown(data.left[1:left.extra])
    data_right<-data.right[1:right.extra]
    data_new1<-c(data_left,data,data_right)
   
    left.zero = rep(0,left.extra)
    right.zero = rep(0,right.extra)
    data_new2<-c(left.zero,data,right.zero)
    

    if      (method == "repeat") data_new<-data_new1
    else if (method == "zero")   data_new<-data_new2
    
    ff<-seq(f1,f2,by=delta.f)
    tt<-seq(t1,t2,by=delta.t)

    length.f<-length(ff)
    length.t<-length(tt)  

    kzft<-array(NA,dim=c(length.t,length.f))
    kzw <-array(NA,dim=c(length.t,length.f))
    em  <-array(NA,dim=c(length.t,length.f))
    pg  <-array(NA,dim=c(length.t,length.f))

    for ( f in (1:length.f) )  { 

    m = round(n/ff[f])
    coef<-coeff.kzft(m,k)
    M<-(m-1)*k+1
    left.half.M<-floor((M-1)/2)
    right.half.M<-ceiling((M-1)/2)
    s<-0:(M-1)
    coefft<-coef*exp(-1i*s*2*pi*ff[f]) 

    for( t in (1:length.t) ) {
 
       kzft[t,f]<-data_new[(tt[t]+left.extra-left.half.M):(tt[t]+left.extra+right.half.M)]%*%coefft
       kzw[t,f]<-2*Re(kzft[t,f])
       em[t,f]<-kzw[t,f]^2
       pg[t,f]<-((abs(kzft[t,f]))^2)*M

    }
}

    lst<-list(kzw=kzw, em=em, pg=pg)
    return(lst)
}


kzww<-function(x,f1=1/length(x),f2=0.5,delta.f=1/length(x),t1=1,t2=length(x),delta.t=1, m, k=3, method="zero"){

updown<-function(x){
    n<-length(x)
    y<-rep(NA,n)

    for ( i in (1:n)) {
     y[i] = x[n-i+1]
    }
    return(y)
}

    data<-as.vector(x)

    left.extra<-floor((m-1)*k/2)
    right.extra<-ceiling((m-1)*k/2)

    left.n.extra<-floor(left.extra/length(data)/2)+1
    right.n.extra<-floor(right.extra/length(data)/2)+1    

    data.left<-rep(c(data,updown(data)),left.n.extra)
    data.right<-rep(c(updown(data),data),right.n.extra)
    
    data_left <-updown(data.left[1:left.extra])
    data_right<-data.right[1:right.extra]
    data_new1<-c(data_left,data,data_right)
   
    left.zero = rep(0,left.extra)
    right.zero = rep(0,right.extra)
    data_new2<-c(left.zero,data,right.zero)
    

    if      (method == "repeat") data_new<-data_new1
    else if (method == "zero")   data_new<-data_new2
    
    ff<-seq(f1,f2,by=delta.f)
    tt<-seq(t1,t2,by=delta.t)

    length.f<-length(ff)
    length.t<-length(tt)  

    kzft<-array(NA,dim=c(length.t,length.f))
    kzw <-array(NA,dim=c(length.t,length.f))
    em  <-array(NA,dim=c(length.t,length.f))
    pg  <-array(NA,dim=c(length.t,length.f))

    coef<-coeff.kzft(m,k)
    M<-(m-1)*k+1
    left.half.M<-floor((M-1)/2)
    right.half.M<-ceiling((M-1)/2)
    s<-0:(M-1)
    omega<-2*pi*ff
    coefft<-coef*exp(-1i*s%o%omega)

    for( t in (1:length.t) ) {
       kzft[t,]<-data_new[(tt[t]+left.extra-left.half.M):(tt[t]+left.extra+right.half.M)]%*%coefft
       kzw[t,]<-2*Re(kzft[t,])
       em[t,]<-kzw[t,]^2
       pg[t,]<-((abs(kzft[t,]))^2)*M

    }

    lst<-list(kzw=kzw, em=em, pg=pg)
    return(lst) 
}

