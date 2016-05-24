
BPSpriorElicit<-function(r0=1,H=1,T00=1,ord=4,G=30,c=0.9){
    hyp<-list(r0=r0,H=H,T00=T00,ord=ord,G=G,c=c)
    # determine spline knots
    hyp$knots<-seq(0,hyp$T00,len=hyp$G)
    extraknots<-cumsum(rep(hyp$knots[2]-hyp$knots[1],hyp$ord-1))
    if(hyp$ord>1) hyp$knots<-c(-extraknots,hyp$knots,hyp$knots[hyp$G]+extraknots)
    # determine spline coefficient variance
    hyp$w<-hyp$H^2
    # determine spline coefficient mean
    hyp$m<-log(hyp$r0)-0.5*hyp$w
    # return hyperparameters
    return(hyp)
} # end BPSpriorElicit
