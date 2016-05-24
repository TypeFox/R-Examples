val4symb<-function(x, FUN=mean, col=c("blue","red"),...){
    FUN <- match.fun(FUN)
    signe<-x-FUN(x,...)
    diam<-abs(x-FUN(x,...))
    color<-ifelse(signe<0,col[1],col[2])
    list(size=diam,col=color)
}
