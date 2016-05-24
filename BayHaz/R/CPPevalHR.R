
CPPevalHR<-function(time,sample){
    ss<-nrow(sample$xi0) # sample size
    hr<-matrix(NA,ss,length(time))
    if (ss>0) for(it in 1:ss){
        # contribution in the origin
        hr[it,]<-sample$xi0[it,1]*sample$hyp$q*pnorm(time,sd=sample$hyp$sd,lower.tail=FALSE)
        # other contributions
        for(j in 1:sample$hyp$F) hr[it,]<-hr[it,]+sample$csi[it,j]*dnorm(time-sample$sgm[it,j],sd=sample$hyp$sd)
    } # end if (ss>0) for(it in 1:ss)
    return(hr) # return hazard rates
} # end CPPevalHR
