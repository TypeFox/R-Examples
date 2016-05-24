Sn.alfa<-function(alfa,n,d,model=c("parametric","semiparametric"), tol = .Machine$double.eps^0.25, maxiter = 1000)
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
    CV<-CV^2
    return(CV)
}

CV.Epidemic.Vn<-function(alfa,d,tol=1.0e-10)
{
    CV<-uniroot(function(x,alfa,d,toler) p.Epidemic.Vn(x,d,tol)-alfa,c(.1,5),alfa=alfa,d=d,toler=tol,tol=0.000001)$root
    return(CV)
}


CV.Epidemic.Wn<-function(alfa,tol=1.0e-7)
{
    CV<-uniroot(function(x,alfa,toler) p.Epidemic.Wn(x,tol)-alfa,c(.1,5),alfa=alfa,toler=tol,tol=0.000001)$root
    return(CV)
}
