ar.bugs <-
function(y, ar.order=1, k=NULL, sim=FALSE, 
                  mean.centre=FALSE, beg=ar.order+1,
                  mean.prior=ar.prior, ar.prior="dnorm(0,1)", tol.prior="dgamma(0.000001,0.000001)", var.prior=NULL, sd.prior=NULL,
                  space=FALSE){
  y<-c(y)
  n<-length(y)
  if(!is.null(k)){
    y<-c(y,rep(NA,k))
  }
  k<-length(y)-max(which(!is.na(y)))
  if(beg<ar.order)
    stop("The value of beg must be at least 1 greater than the number of lags")
  if(!is.null(var.prior) | !is.null(sd.prior)){
    tol.prior<-NULL
  }
  if(length(c(tol.prior,var.prior,sd.prior))>1)
    stop("Only one of tol.prior, var.prior or sd.prior should be given. Set others to null")
  
  bug<-c("model{","")
  #likelihood
  lik<-c("#likelihood",
         paste0("for(t in ",beg,":",n+k,"){"),
         "\ty[t] ~ dnorm(y.mean[t], isigma2)",
         "}")
  bug<-c(bug, lik)
  #ymean
  ymean<-c("#mean",
           paste0("for(t in ",beg,":",n+k,"){"),
           y.mean<-c("\ty.mean[t] <- 0",
                     "}")
  )
  if(ar.order==0 & mean.centre==T)  ymean[3]<-"\ty.mean[t] <- phi0"
  if(ar.order!=0 & mean.centre==F)  ymean[3]<-paste0("\ty.mean[t] <- ",paste0("phi",1:ar.order,"*y[t-",1:ar.order,"]",collapse=" + "))
  if(ar.order!=0 & mean.centre==T)  ymean[3]<-paste0("\ty.mean[t] <- phi0 + ",paste0("phi",1:ar.order,"*(y[t-",1:ar.order,"]-phi0)",collapse=" + "))
  bug<-c(bug, ymean, "")
  
  #prior
  if(!is.null(tol.prior)){
    ar.priors<-c(paste0("phi",0:ar.order," ~ ",ar.prior), 
                 paste0("isigma2 ~ ",tol.prior),
                 "sigma <- pow(isigma2,-0.5)")
  }
  if(!is.null(var.prior)){
    ar.priors<-c(paste0("phi",0:ar.order," ~ ",ar.prior), 
                 paste0("sigma2 ~ ",var.prior),
                 "isigma2 <- pow(sigma2,-1)")
  }
  if(!is.null(sd.prior)){
    ar.priors<-c(paste0("phi",0:ar.order," ~ ",ar.prior), 
                 paste0("sigma ~ ",sd.prior),
                 "isigma2 <- pow(sigma,-2)")
  }
  ar.priors<-ar.priors[-1]
  if(mean.centre==T)  ar.priors<-c(paste0("phi0 ~ ",mean.prior),ar.priors)
  bug<-c(bug,"#priors",ar.priors,"")
  
  #forecast
  forc<-NULL
  if(k!=0){
    forc<-c("#forecast",
            paste("for(t in ",n+1,":",n+k,"){",sep=""),
            "\ty.new[t] <- y[t]",
            "}",
            "")
    bug<-c(bug,forc)
  }  
  
  #simulation
  ysim<-NULL
  if(sim==TRUE){
    ysim<-c("#simulation",
            "isigma2.c <- cut(isigma2)",
            paste("for(t in ",beg,":",n,"){",sep=""),
            "\ty.mean.c[t] <- cut(y.mean[t])",
            "\ty.sim[t] ~ dnorm(y.mean.c[t],isigma2.c)",
            "}",
            "")
    bug<-c(bug,ysim)
  }
  bug<-c(bug,"}","")
  #print.tsbugs(list(bug=bug))
  
  if(space==FALSE){
    bug<-bug[-nchar(bug)!=0]
    if(length(grep("#mean", bug))>0)
      bug<-bug[-grep("#mean", bug)]
  }
    
  p1<-grep("#likelihood",bug)
  p2<-grep("#priors",bug)
  if(k!=0 & sim==TRUE){
    p3<-grep("#forecast",bug); p4<-grep("#simulation",bug)
  }
  if(k!=0 & sim==FALSE){
    p3<-grep("#forecast",bug); p4<-length(bug)
  }
  if(k==0 & sim==TRUE){
    p3<-grep("#simulation",bug); p4<-p3
  } 
  if(k==0 & sim==FALSE){
    p3<-length(bug); p4<-p3
  } 
  p5<-length(bug)
  
  bug<-list(bug=bug,
            data=list(y=y),
            info=list(n=n,k=k,nh=n+k,beg=beg,
                      args=mget(names(formals()),sys.frame(sys.nframe()))[-1],
                      variance="CV",
                      likelihood=p1:(p2-1),
                      priors=p2:(p3-1),
                      forecast=NULL,
                      simulation=NULL))
  if(p3!=p4)  bug$info$forecast<-p3:(p4-1)
  if(p4!=p5)  bug$info$simulation<-p4:(p5-1)
  class(bug)<-"tsbugs"
  return(bug)
}
