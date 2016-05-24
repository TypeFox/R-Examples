
BPSevalHR<-function(time,sample){
    ss<-nrow(sample$eta) # sample size
    hr<-matrix(NA,ss,length(time))
    B<-splineDesign(knots=sample$hyp$knots,x=time,ord=sample$hyp$ord)
    if (ss>0) for(it in 1:ss) hr[it,]<-exp(B%*%sample$eta[it,])
    return(hr) # return hazard rates
} # end BPSevalHR
