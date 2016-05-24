
CPPpriorSample<-function(ss=1,hyp=CPPpriorElicit()){
    sampri<-list(hyp=hyp,sgm=matrix(NA,ss,hyp$F),xi0=matrix(NA,ss,1),csi=matrix(NA,ss,hyp$F))
    if (ss>0) for(it in 1:ss){
        #   hazard sources (jump-times)
        sampri$sgm[it,]<-cumsum(rexp(hyp$F,rate=hyp$q))
        #   hazard contribution in the origin (jump-size in the origin)
        sampri$xi0[it,1]<-rgamma(1,shape=hyp$a,rate=hyp$b)
        #   other hazard contributions (other jump-sizes)
        sampri$csi[it,]<-rgamma(hyp$F,shape=hyp$a,rate=hyp$b)
    } # end if (ss>0) for(it in 1:ss)
    return(sampri)
} # end CPPpriorSample
