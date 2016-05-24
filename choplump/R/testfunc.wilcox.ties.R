`testfunc.wilcox.ties` <-
function(choplist,RM,NTIES,M){
    ZM<-choplist$ZM
    a<-choplist$a
    b<-choplist$b
    n0<-length(ZM[ZM==0])+a
    n1<-M+a+b-n0
    r0<- (a+b)-(a+b-1)/2
    RM<-RM+a+b
    NTIES<-c(a+b,NTIES)
    STATISTIC <- sum(RM[ZM==0])+a*r0 - n0 * (n0 + 1)/2
    SIGMA <- sqrt((n0 * n1/12) * ((n0 + n1 + 1) - 
            sum(NTIES^3 - NTIES)/((n0 + n1) * (n0 + n1 - 
              1))))
    out<-  (STATISTIC - n0 * n1/2)/SIGMA
    return(out) 
}

