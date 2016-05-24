boldContrast <-
function(B0,te,plot=TRUE,random=TRUE,alpha=.38,hct=.4,cbf.base=55,chi=1.8e-7,e.base=.4,w0=2.675e8,e.act=.1,cbf.act=2*cbf.base){
  M0<-100
  cbv.base <- 0.8 * cbf.base^alpha
  vol.base <- cbv.base/100
  cbv.act<- 0.8 * cbf.act^alpha
  vol.act<-cbv.act/100
  r2.si <- (1.74 * B0 + 7.77)
  chi.si<-chi*4*pi
  w0.si = w0/(2*pi)
  ######################################
  ######## CALCULATE BASE T2* ##########
  ######################################
  if(random==FALSE){
    f.shift.si<-w0.si * B0 * chi.si * hct * (e.base) * 2 * pi}else{
      f.shift.si <- w0.si * B0 * chi.si * hct * (e.base) * 4 * pi/3
    }  
  
  r.si <- vol.base * f.shift.si
  
  r2.star.base <- r.si + r2.si
  t2.star.base.si <- 1/r2.star.base * 1000
  #######################################
  ######## CALCULATE ACTIVE T2* #########
  #######################################
  if(random==FALSE){
    r2.star.act <- r2.si+(vol.act * w0.si * B0 * chi.si * hct * (e.act) * 2 * pi)}else{
      r2.star.act <- r2.si+(vol.act * w0.si * B0 * chi.si * hct * (e.act) * 4 * pi/3)
    }
  t2.star.act.si<-1/r2.star.act*1000
  #######################################
  ########## SIGNAL CONTRAST ############
  #######################################
  
  t<-0:250
  signal.t2.star.base<- M0*(exp(-t/t2.star.base.si))  ## Signal displays an exponential dependence on TE relative to base T2* 
  signal.t2.star.act<- M0*(exp(-t/t2.star.act.si))    ## Signal displays an exponential dependence on TE relative to active T2*
  sig.dif<-signal.t2.star.act-signal.t2.star.base       ## Subtract active signal from base signal to measure bold sensitivity 
  bs<-sig.dif/max(sig.dif)*100                          ## convert BOLD sensitivity to percentage scale
  study.bs<-(sig.dif[te+1]/max(sig.dif))*100 
  
  if(plot==TRUE){  
    signal.t2.star.base<- M0*(exp(-t/t2.star.base.si))  ## Signal displays an exponential dependence on TE relative to base T2* 
    signal.t2.star.act<- M0*(exp(-t/t2.star.act.si))    ## Signal displays an exponential dependence on TE relative to active T2*
    sig.dif<-signal.t2.star.act-signal.t2.star.base        ## Caclculate the expected BOLD contrast accross TEs
    bs<-sig.dif/max(sig.dif)*100                           ## Convert to sensativity(%)
    title<-paste("BOLD Signal as a Function of T2* and TE (B0 = ", 
                 as.character(B0),"T)", sep=" ")           ## Create title of graph
    
    plot(t,signal.t2.star.base,lwd=7,                      ## Plot it 
         type="l",col="red",ylim=c(0,100),
         xlab="TE (ms)",ylab="Signal (%)",
         main=title)
    
                        ## Plot BOLD Sensitivity
    
    x<-c(0,t[length(t)])                                  ## Do the Nice Lines
    y<-matrix(rep(seq(0,100,100/5),each=2),
              ncol=2,byrow=TRUE)
    for(i in 1:nrow(y)){
      lines(x,y[i,],col="grey")
    }            
    lines(t,signal.t2.star.act,lwd=7,type="l",col="blue") ## Plot Activation T2*
    lines(t,bs,col="green",lwd=7)   
    lines(t,signal.t2.star.base,lwd=7,type="l",col="red") 
    legend(x="topright",                                  ## Create the Legend
           legend=c("T2* Activation",
                    "T2* Baseline",
                    "BOLD Sensitivity(%)"),
           col=c("blue","red","green"),lwd=2,bg="white")
    
  }
  return(sig.dif[te+1])
}
