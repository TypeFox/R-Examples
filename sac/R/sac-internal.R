dvalue.epidemic<-function(x,tol=1.0e-7)
{
    delta<-1
    d.value<-0; k<-1
    while(abs(delta)>tol){
        delta<-(4*k^2*x)*(4*k^2*x^2-3)*exp(-2*k^2*x^2)
        d.value<-d.value+delta
        k<-k+1
    }
    d.value<-2*d.value
    return(d.value)
}

Deltan.ModelTest<-function(x,k,m,Alpha,Beta)
{
    if(is.vector(x)){
        n<-length(x)
        z<-sort(x)
        x<-matrix(x,n,1)
    }
    if(!is.matrix(x)) x<-as.matrix(x)
    if(is.matrix(x)){
        n<-nrow(x)
        z<-sort(x[,1])
    }
    p<-NULL
    eF<-NULL; sF<-NULL; eG<-NULL; sG<-NULL
    eF[1]<-0; sF[1]<-0; eG[1]<-0; sG[1]<-0
    p<-1/((m-k)*exp(Alpha+x%*%Beta)+(n-m+k))
    for(i in 1:n){
        if(m<n) xx<-c(x[1:k,1],x[(m+1):n,1])
        else xx<-x[1:k,1]
        eF[i+1]<-sum(xx<=z[i])/(n-m+k)
        eG[i+1]<-sum(x[(k+1):m,]<=z[i])/(m-k)
        sF[i+1]<-sum(p*(x[,1]<=z[i]))
        sG[i+1]<-sum(exp(Alpha+x%*%Beta)*p*(x[,1]<=z[i]))
    }
    Delta<-max((k+n-m)*abs(sF-eF)+(m-k)*abs(sG-eG))/sqrt(n)
    return(Delta)
}
G.alfa<-function(alfa,x,n,d,model="semiparametric")
{
    if(model=="semiparametric"){
        h<-(log(n)/n)^((1-alfa)*(d+1/2)); l<-h  # semiparametric
    }
    else{
        h<-(log(n))^((1-alfa)*(1+1/2))/n^1; l<-h  # parametric
    }
    return( dchisq(x^2,d)*((x^2-d)*log((1-l)*(1-h)/(l*h))+4)-alfa)
}

Z.alfa<-function(alfa,n,d)
{
    t.alfa<--log(log(1/sqrt(1-alfa)))
    x<-log(n)
    Dd<-2.*log(x)+0.5*d*log(log(x))-lgamma(0.5*d)
    return((t.alfa+Dd)/sqrt(2*log(log(n))))
}

Zn.alfa<-function(alfa,n,d,model=c("parametric","semiparametric"), tol = 1.0e-7, maxiter = 1000)
{
    model <- match.arg(model)
    if(model=="parametric") h<-(log(n))^((1-alfa)*(1+1/2))/n^1  # parametric
    else h<-(log(n)/n)^((1-alfa)*(d+1/2)) # semiparametric
    tempfun<-function(x,alfa,d,h) dchisq(x^2,d)*((x^2-d)*log((1-h)^2/(h*h))+4)-alfa
    d.tempfun<-function(x,d,h)
    {
        A<-log((1-h)^2/(h*h))
        dchisq(x^2,d)*(-8-A*d^2-4*x^2-A*x^4+2*d*(2+A*(1+x^2)))/x
    }
    CV<-Z.alfa(alfa,n,d)
    del<-abs(tempfun(CV,alfa,d,h))
    it<-1
    while(it<=50 || del>tol)
    {
        CV.nu<-CV-tempfun(CV,alfa,d,h)/d.tempfun(CV,d,h)
        del<-abs(tempfun(CV.nu,alfa,d,h))
        CV<-CV.nu
        it<-it+1
    }
    return(CV)
}
