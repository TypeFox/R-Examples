`testfunc.wilcox` <-
function(ZM,n1,n0,RM){
    a<-n0- length(ZM[ZM==0])
    b<-n1-length(ZM[ZM==1])
    r0<- (a+b)-(a+b-1)/2
    RM<-RM+a+b
    out <- sum(RM[ZM==0])+a*r0 
    return(out) 
}

