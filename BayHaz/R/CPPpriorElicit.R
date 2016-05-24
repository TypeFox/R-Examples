
CPPpriorElicit<-function(r0=1,H=1,T00=1,M00=1,extra=0){
    hyp<-list(r0=r0,H=H,T00=T00,M00=M00)
    # select exponential compound distribution
    hyp$a<-1
    # determine bandwidth (Gaussian kernel standard deviation)
    hyp$sd<-sqrt((hyp$T00)^2/(8*hyp$M00^2))
    # determine underlying Poisson intensity (hazard sources per time unit)
    hyp$q<-(2^(1/2)*hyp$M00*(1+hyp$a^(-1)))/(pi^(1/2)*hyp$H^2*hyp$T00)
    # determine compound distribution rate parameter (events per time unit)
    hyp$b<-hyp$q*hyp$a/hyp$r0
    # determine "a priori" maximum useful number of hazard sources (not counting the one in the origin)
    hyp$F<-qpois(0.95,hyp$q*(hyp$T00+3*hyp$sd))+extra
    # return hyperparameters
    return(hyp)
} # end CPPpriorElicit
