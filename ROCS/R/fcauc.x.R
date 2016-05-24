fcauc.x <-
function(x0, x1, FDR.cut=0.2)
{
    if(median(x1)<median(x0))
    {
        x1<- -x1
        x0<- -x0
    }
    
    x<-c(x0, x1)
    y<-c(rep(0, length(x0)), rep(1,length(x1)))
    
    o<-order(x)
    y<-y[o]
    x<-x[o]
    
    l0<-length(x0)
    l1<-length(x1)
    l<-l0+l1
    
    bb<-cumsum(y[1:(l-1)])	#count of true positive being claimed negative
    aa<-1:(l-1) - bb	#count of true negative claimed negative
    dd<-l1 - bb	#count of true positive claimed positive
    cc<-l - 1:(l-1) - dd	#count of true negative claimed positive
    
    FP<-cc/l0
    TP<-dd/l1
    FDR<-cc/(cc+dd)
    FDR<- -cummax(-FDR)
    TDR<-1-FDR
    
    fcauc<-fcauc.fptp(FP, TP, TDR, FDR.cut=FDR.cut)
    fcauc
}
