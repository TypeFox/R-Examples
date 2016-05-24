TAR.thres<-function(ay,p1,p2,ph.1,ph.2,sig.1,sig.2,lagd,thres,step.r=0.02,bound,lagp1,lagp2,constant=1,thresVar){   ## step.r is the step size of the MH sampling.
                                                                                     ## bound is the hyper-parameter of a Gamma prior
new.r<- thres+step.r*rnorm(1,mean=0,sd=1)    ## Sampling a candidate threshold value

repeat{                                         ## Check whether the new threshold value is located on U(a,b).
if((new.r< bound[1])|(new.r> bound[2])){        ## If not, repeat to sample a new value.
new.r<- thres+step.r*rnorm(1,mean=0,sd=1)}
else break
}

## The random walk M-H algorithm
if (!missing(thresVar)){
old.lik<- TAR.lik(ay,p1,p2,ph.1,ph.2,sig.1,sig.2,lagd,thres,lagp1,lagp2,constant=constant,thresVar)
new.lik<- TAR.lik(ay,p1,p2,ph.1,ph.2,sig.1,sig.2,lagd,new.r,lagp1,lagp2,constant=constant,thresVar)
}
else{
old.lik<- TAR.lik(ay,p1,p2,ph.1,ph.2,sig.1,sig.2,lagd,thres,lagp1,lagp2,constant=constant)
new.lik<- TAR.lik(ay,p1,p2,ph.1,ph.2,sig.1,sig.2,lagd,new.r,lagp1,lagp2,constant=constant)
}
if((new.lik-old.lik)>log(runif(1))){r.count=1}   ## To determine whether update a new value or not.
else{ new.r<- thres; r.count<- 0 }

return(c(r.count,new.r))
}

