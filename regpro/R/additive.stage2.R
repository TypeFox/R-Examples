additive.stage2<-function(arg,x,y,h=1,kernel="gauss",B=2)
{
d<-length(arg)

if (kernel=="gauss") ker<-function(t){ return( exp(-t^2/2) ) }
if (kernel=="uniform") ker<-function(t){ return( (abs(t) <= 1) ) }

ynow<-y
resu<-0
ssr<-matrix(0,d,1)
for (ii in 1:B){
    for (jj in 1:d){
        xjj<-x[,jj]
        xjj<-matrix(xjj,length(xjj),1)
        argu<-matrix(arg[jj],length(xjj),1)
        w<-ker((xjj-argu)/h)/h^d
        w<-w/sum(w)
        yhat<-w%*%ynow
        ssr[jj]<-sum((yhat-ynow)^2)
    }  
    dstar<-which.min(ssr)
    arg.now<-arg[dstar]
    x.now<-x[,dstar]

    w<-kernesti.weights(arg.now,x.now,h=h,kernel=kernel)
    w<-matrix(w,length(w),1)    
    ynow<-matrix(ynow,length(ynow),1)
    notna<-(!is.na(ynow))
    w<-notna*w
    w<-w/sum(w)
    mu<-w*ynow
    neweva<-sum(mu,na.rm=TRUE)

    resu<-resu+neweva    
    residu<-1 
    ynow<-residu        
}

return(resu)
}







