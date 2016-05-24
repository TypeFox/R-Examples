p.OneChange<-function(n,d,Sn)
{
    A<-sqrt(2*log(log(n)))
    Dd<-2*log(log(n))+0.5*d*log(log(log(n)))-lgamma(d/2)
    p.value<-1-exp(-2*exp(-(A*sqrt(Sn)-Dd)))
    
    tempfun<-function(alfa,x,n,d){
        h<-(log(n)/n)^((1-alfa)*(d+1/2)) # semiparametric
        dchisq(x,d)*((x-d)*log((1-h)^2/(h*h))+4)-alfa
    }
    p.value<-uniroot(tempfun,c(0,p.value),x=Sn,n=n,d=d,tol=1e-7)$root
    return(p.value)
}


p.Epidemic.Vn<-function(Vn,d,tol=1.0e-10)
{
    delta<-1; p.value<-0; 
    if(d==1){
        i<-1; 
        while(abs(delta)>tol){
            delta<-2*(-1)^(i+1)*exp(-pi^2*i^2*Vn)
            p.value<-p.value+delta
            i<-i+1
        }
    }
    if(d==2){
        k<-1;
        p1<-0; p2<-0
        del1<-1;del2<-1
        while(abs(delta)>tol){
            del1<-(pi^2*k^2*Vn+1)*exp(-pi^2*k^2*Vn)
            p1<-p1+del1
            for(i in 1:k){
                j<-k+1
                del2<-(-1)^(i+j)*(j^2*exp(-pi^2*i^2*Vn)-i^2*exp(-pi^2*j^2*Vn))/(j^2-i^2)
                p2<-p2+del2
            }
            delta<-4*(del1+del2)
            p.value<-4*(p1+p2)
            k<-k+1
        }
    }
    return(p.value)
}


p.Epidemic.Wn<-function(Wn,tol=1e-7)
{
    delta<-1
    p.value<-0; k<-1
    while(abs(delta)>tol){
        delta<-(4*k^2*Wn-1)*exp(-2*k^2*Wn)
        p.value<-p.value+delta
        k<-k+1
    }
    p.value<-2*p.value
    return(p.value)
}
