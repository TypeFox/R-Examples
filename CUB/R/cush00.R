# @title CUSH model without covariates
# @description Estimate and validate a CUSH model for given ordinal responses, without covariates.
# @usage cush00(m, ordinal, shelter, makeplot)
# @param m Number of ordinal categories
# @param ordinal Vector of ordinal responses
# @param shelter Category corresponding to the shelter choice
# @param makeplot Logical: if TRUE, the algorithm returns a graphical plot comparing fitted 
# probabilities and observed relative frequencies 
# @aliases cush00
# @return An object of the class "CUSH"
# @import stats graphics
#' @keywords internal

cush00 <-
function(m,ordinal,shelter,makeplot){
  tt0<-proc.time()
  freq<-tabulate(ordinal,nbins=m); n<-length(ordinal); aver<-mean(ordinal);
  fc<-freq[shelter]/n
  deltaest<-max(0.01,(m*fc-1)/(m-1))   ### sufficient unbiased estimator
  esdelta<-sqrt((1-deltaest)*(1+(m-1)*deltaest)/(n*(m-1)))
  varmat<-esdelta^2
  wald<-deltaest/esdelta
  loglik<-loglikcush00(m,ordinal,deltaest,shelter)
  AICCUSH<- -2*loglik+2
  BICCUSH<- -2*loglik+log(n)
  llunif<- -n*log(m); csisb<-(m-aver)/(m-1);
  llsb<-loglikcub00(m,freq,1,csisb)
  nonzero<-which(freq!=0)
  logsat<- -n*log(n)+sum((freq[nonzero])*log(freq[nonzero]))
  devian<-2*(logsat-loglik)
  LRT<-2*(loglik-llunif)
  ###########
  cat("==================================================================","\n")
  cat(">>ML estimate of a CUSH model (September 2015) <<","\n") 
  cat("==================================================================","\n")
  cat("n = ", n,"    m =",m,"   shelter =",shelter,"\n")
  cat("==================================================================","\n")
  cat("ML estimate of delta  =",round(deltaest,digits=5),"\n")
  cat("Standard error        =",round(esdelta,digits=5),"\n")
  cat("Estimate/Stand.err.   =",round(wald,digits=5),  "\n")
  cat("LRT test              =",round(LRT,digits=5),    "\n")
  cat("==================================================================","\n")
  cat("Log-lik(delta^)    =",round(loglik,digits=8),"\n")
  cat("Log-lik(saturated) =",round(logsat,digits=8),"\n")
  cat("Deviance           =",round(devian,digits=8),"\n")
  cat("------------------------------------------------------------------","\n")
  cat("Log-lik(UNIFORM)         =",round(llunif,digits=8),"\n")
  cat("Log-lik(Shifted-BINOMIAL)=",round(llsb,digits=8),"\n")
  cat("------------------------------------------------------------------","\n")
  cat("AIC-CUSH   =",round(AICCUSH,digits=8),"\n")
  cat("BIC-CUSH   =",round(BICCUSH,digits=8),"\n")
  cat("==================================================================","\n")
  ### Observed relative frequencies and estimated probabilities
  theorpr<-deltaest*ifelse(seq(1,m)==shelter,1,0)+(1-deltaest)/m
  pearson<-((freq-n*theorpr))/sqrt(n*theorpr)
  X2<-sum(pearson^2)
  relares<-(freq/n-theorpr)/theorpr
  diss00<-dissim(theorpr,freq/n)
  FF2<-1-diss00
  LL2<-1/(1+mean((freq/(n*theorpr)-1)^2))
  II2<-(loglik-llunif)/(logsat-llunif)
  cat("Pearson Fitting measure    ==>  X^2 =",X2,"(p-val.=",1-pchisq(X2,m-3),")","\n")
  cat("Lik-based fitting measure  ==>  L^2 =",LL2,"\n")
  cat("Relative Log-lik index     ==>  I   =",round(II2,digits=5),"\n")
  cat("F^2 fitting measure        ==>  F^2 =",round(FF2,digits=5),"\n")
  cat("Normed Dissimilarity index ==>  Diss=",diss00,"\n")
  cat("==================================================================","\n")
  cat("(R=r) Observed CUB-prob","Pearson","Relative res.","\n")
  stampa<-cbind(1:m,freq/n,theorpr,pearson,relares)
  print(stampa,digits=5)
  cat("==================================================================","\n")
  ### Check on plotting
  ### if(makeplot==TRUE) DO PLOT; else (makeplot==FALSE) ===> NOPLOT
  ### ===>>> ## stringtitle=match.call()[[2]]; ### it writes "ordinal"
  if(makeplot==TRUE){
    ################### GRAFICI sovrapposti #########
    par(mar=c(4,4,2.5,1)+0.1) ;               ### reset standard margins
    par(mfrow=c(2,1));            ### ripristina l'area del grafico
    ###### Distributions
    stringtitle="CUSH model"; 
    plot(cbind(1:m,1:m),cbind(theorpr,(freq/n)),las=1,
         main=paste(stringtitle,  "     (Diss =",round(diss00,digits=4),")"),
         xlim=c(1,m),ylim=c(0.0,1.1*max(theorpr,(freq/n))),
         xlab="", ylab=""); #xlab=Ordinal values of R=1,2,...,m
    points(1:m,theorpr,pch=21,cex=2,lwd=2.0,type="b",lty=3,col="red");
    points(1:m,freq/n,pch=19,cex=1.5,lwd=1.5);
    abline(h=0);
    ###### Log-likelihood
    vettdelta<-seq(0,deltaest*1.1,by=0.001)  #vettdelta=seq(0,0.3,by=0.001)
    vettelle<-n*((1-fc)*log(1-vettdelta)+fc*log(1+(m-1)*vettdelta)-log(m))
    minelle<-min(vettelle)
    estremo<-loglik-2
    dove<-(vettelle>estremo)
    interv<-vettdelta[dove==TRUE]
    deltamin<-interv[1];deltamax<-interv[length(interv)]
    titolo<-paste("Log-likelihood function", "             delta-estimate = ", round(deltaest,digits=3))
    plot(vettdelta,vettelle,main=titolo,ylab="",xlab=expression(delta),ylim=c(llunif,logsat),
         pch=19,cex.main=0.9,lwd=1,col="blue")
    ### arrows(deltamin,minelle,deltamax,minelle,code=3,angle=90,length=0.10,lty=3,lwd=1.2,col="red")
    points(c(0,0,0),c(llunif,loglik,logsat),pch=19,cex=1,col="red")
    text(0,llunif,labels="Uniform",pos=4,offset=0.5,font=4,cex=0.7)
    text(0,loglik,labels="CUSH",pos=4,offset=0.5,font=4,cex=0.7)
    text(0,logsat,labels="Saturated",pos=4,offset=0.5,font=4,cex=0.7)
    text(deltamin,minelle,labels=as.character(deltamin),font=4,cex=0.8)
    text(deltamax,minelle,labels=as.character(deltamax),font=4,cex=0.8)
    ### abline(h=estremo,lwd=0.8,col="red")
    points(deltaest,minelle,col="green3",pch=8,cex=1.5,lwd=3)
    par(mar=c(5,4,4,2)+0.1)                ### reset standard margins
    par(mfrow=c(1,1))                      ### ripristina l'area del grafico
    ####################################################################
  }
  #####################################################################
  # Assignments as global variables: assign('name',value,pos=1)
  #####################################################################
  #   assign('deltaest',deltaest,pos=1)
  #   assign('loglik',loglik,pos=1)
  #   assign('varmat',varmat,pos=1)
  ####################################
  durata<-proc.time()-tt0;durata<-durata[1];
  cat("=======================================================================","\n")
  cat("Elapsed time=",durata,"seconds","=====>>>",date(),"\n")
  cat("=======================================================================","\n","\n")
  results<-list('estimates'=round(deltaest,digits=5), 'loglik'=loglik,'varmat'=varmat,'BIC'= round(BICCUSH,digits=8))
}
