`saws` <-
function(x,test=diag(p),beta0=matrix(0,p,1),conf.level=0.95,method=c("d3","d5","d1","d2","d4","dm"),
    bound=.75){
    p<-length(x$coefficients)
    out<-sawsCalc(x$coefficients,x$u,x$omega,test,beta0,conf.level,method,bound)
    if (!is.null(x$originalCall)) out$originalCall<-x$originalCall
    else if (!is.null(x$call)) out$originalCall<-x$call
    out
}
