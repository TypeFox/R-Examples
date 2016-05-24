
BPSpriorSample<-function(ss=1,hyp=BPSpriorElicit()){
    ns<-length(hyp$knots)-hyp$ord # number of splines
    sampri<-list(hyp=hyp,eta=matrix(NA,ss,ns))
    sampri$eta[,1]<-rnorm(ss,mean=hyp$m,sd=sqrt(hyp$w))
    for(g in seq(1,ns-1))
        sampri$eta[,g+1]<-hyp$c*sampri$eta[,g]+rnorm(ss,mean=(1-hyp$c)*hyp$m,sd=sqrt((1-hyp$c^2)*hyp$w))
    return(sampri)
} # end BPSpriorSample
