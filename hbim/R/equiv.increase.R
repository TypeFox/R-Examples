`equiv.increase` <-
function(x1,y1,x2,y2,e1,xlog=TRUE){
    a1<-approx(y1,x1,xout=e1)$y
    e2<-approx(x2,y2,xout=a1)$y
    a2<-approx(y1,x1,xout=e2)$y
    if (xlog) equiv.increase<- 10^(a2-a1)
    else equiv.increase<- a2/a1
    out<-list(a1=a1,e2=e2,a2=a2,e1=e1,equiv.increase=equiv.increase)
    return(out)
}

