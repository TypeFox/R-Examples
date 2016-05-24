`calc.foldrange` <-
function(n,lower,upper,conf.level=.95){
    alpha<-1-conf.level
    s.t<- (sqrt(n)*(log10(upper)-log10(lower)))/(2*qt(1-alpha/2,n-1))
    s.z<-(sqrt(n)*(log10(upper)-log10(lower)))/(2*qnorm(1-alpha/2))
    range.t<- 10^(( qnorm(1-alpha/2) - qnorm(alpha/2) )*s.t )
    range.z<- 10^(( qnorm(1-alpha/2) - qnorm(alpha/2) )*s.z )
    col.names<-c("n",paste("lower",100*conf.level,"% cl"),
           paste("upper",100*conf.level,"% cl"),
           "standard deviation (log10 scale) by t",
           "standard deviation (log10 scale) by Z",
           paste(100*conf.level,"% fold-range, s estimated by t"),
           paste(100*conf.level,"% fold-range, s estimated by Z"))
    out<-data.frame(n=n,lower=lower,upper=upper,s.byt=s.t,s.byz=s.z,
        foldrange.byt=range.t,foldrange.byz=range.z)
 
    #out<-list(out,row.names)
    out
}

