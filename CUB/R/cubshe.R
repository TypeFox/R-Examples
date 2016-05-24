# @title Main function for CUB models with a shelter effect
# @description Estimate and validate a CUB model with a shelter effect.
# @aliases cubshe
# @usage cubshe(m, ordinal, shelter, maxiter, toler, makeplot)
# @param m Number of ordinal categories
# @param ordinal Vector of ordinal responses
# @param shelter Category corresponding to the shelter choice
# @param maxiter Maximum number of iterations allowed for running the optimization algorithm 
# @param toler Fixed error tolerance for final estimates 
# @param makeplot Logical: if TRUE (default), the algorithm returns a graphical plot comparing fitted 
# probabilities and observed relative frequencies
# @return An object of the class "CUB"
# @import stats
# @references
# Iannario M. (2012). Modelling \emph{shelter} choices in a class of mixture models for ordinal responses,  
# \emph{Statistical Methods and Applications}, \bold{21}, 1--22
#' @keywords internal 



cubshe <-
function(m,ordinal,shelter,maxiter,toler,makeplot){
  tt0<-proc.time()
  serie<-1:m; freq<-tabulate(ordinal,nbins=m); n<-sum(freq); 
  aver<-mean(ordinal); varcamp<-mean(ordinal^2)-aver^2;
  dd<-ifelse(serie==shelter,1,0)
  #####################################################################
  vett<-inibest(m,freq)
  pai1<-vett[1]; csi<-vett[2];
  fc<-freq[shelter]
  deltaini<-max(0.01,(m*fc-1)/(m-1)) #deltaini=runif(1,0,0.3)  
  pai2<-max(0.01, 1-deltaini-pai1)
  ################################################################
  loglik<-loglikcubshe(m,freq,pai1,pai2,csi,shelter)
  # ********************************************************************
  # ************* E-M algorithm for CUBSHE *****************************
  # ********************************************************************
  nniter<-1
  while(nniter<=maxiter){
    likold<-loglik
    bb<-probbit(m,csi)
    tau1<-pai1*bb
    tau2<-pai2*(1/m)
    denom<-tau1+tau2+(1-pai1-pai2)*dd
    tau1<-tau1/denom
    tau2<-tau2/denom        
    tau3<-1-tau1-tau2
    numaver<-sum(serie*freq*tau1)
    denaver<-sum(freq*tau1)
    averpo<-numaver/denaver
    pai1<-sum(freq*tau1)/n    #updated pai1 estimate
    pai2<-sum(freq*tau2)/n    #updated pai2 estimate
    csi<-(m-averpo)/(m-1)     #updated csi estimate
    
    if(csi<0.001){
      csi<-0.001;nniter<-maxiter-1;
    }
    
    loglik<-loglikcubshe(m,freq,pai1,pai2,csi,shelter)
    liknew<-loglik
    testll<-abs(liknew-likold)
    #print(cbind(nniter,testll,pai1,pai2,csi,loglik));
    if(testll<=toler) break else {loglik<-liknew}
    nniter<-nniter+1
  }
  
  ######
  if(csi>0.999) csi<-0.99   ###????????### to avoid division by 0 !!!
  if(csi<0.001) csi<-0.01  ###????### to avoid division by 0 !!!
  if(pai1<0.001) pai1<-0.01   ###?????         ### to ensure identifiability !!!
  
  ####################################################################
  # Compute asymptotic standard errors of ML estimates
  ####################################################################
  varmat<-varcovcubshe(m,pai1,pai2,csi,shelter,n)
  nomi<-rbind("pai1","pai2","csi")
  stime<-c(pai1,pai2,csi)
  nparam<-length(stime)
  delta<-1-pai1-pai2
  paistar<-pai1/(pai1+pai2)
  if (isTRUE(varmat==matrix(NA,nrow=nparam,ncol=nparam))==TRUE){
    ddd<-cormat<-matrix(NA,nrow=nparam,ncol=nparam)
    ICOMP<-trvarmat<-NA
    esdelta<-pvaldelta<-espaistar<-pvalpaistar<-NA
    errstd<-wald<-pval<-rep(NA,nparam)
  } else {
    ddd<-diag(sqrt(1/diag(varmat)))
    esdelta<-sqrt(varmat[1,1]+varmat[2,2]+2*varmat[1,2])
    pvaldelta<-round(2*(1-pnorm(abs(delta/esdelta))),20)
    espaistar<-sqrt((pai1^2*varmat[2,2]+pai2^2*varmat[1,1]-2*pai1*pai2*varmat[1,2]))/(pai1+pai2)^2
    pvalpaistar<-round(2*(1-pnorm(abs(paistar/espaistar))),20)  
    trvarmat<-sum(diag(varmat))
    ICOMP<- -2*loglik + nparam*log(trvarmat/nparam) - log(det(varmat))
    errstd<-sqrt(diag(varmat))  
    wald<-stime/errstd
    pval<-round(2*(1-pnorm(abs(wald))),20)
  }
  ####################################################################
  # Print CUBSHE results of ML estimation  
  ####################################################################
  cat("\n")
  cat("=======================================================================","\n")
  cat("============== CUB Program: version 4.0 (September 2015) =============","\n")
  cat("=======================================================================","\n")
  cat("==>>> C U B + SHELTER model  <<<=====  ML-estimates via E-M algorithm  ","\n")
  cat("=======================================================================","\n")
  cat(" m=", m," shelter=", shelter," Sample size: n=", n,
      " *** Iterations=",nniter,"Maxiter=",maxiter,"\n")
  cat("=======================================================================","\n")
  cat("parameters  ML-estimates  stand.errors    estimates/stand.errors       ","\n")
  cat("=======================================================================","\n")
  for(i in 1:3){
    cat(nomi[i],"        ",round(stime[i],5),"        ",round(errstd[i],5),"        ",round(wald[i],5),"      ","\n")
  }
  cat("=======================================================================","\n")
  cat("delta","       ",delta, "    ",esdelta,  "       ",delta/esdelta,"     ","\n")
  cat("(1-pai1-pai2)","\n")
  cat("=======================================================================","\n")
  cat("paistar","      ",paistar,"      ",espaistar,"        ",paistar/espaistar,"       ","\n")
  cat("[pai1/(pai1+pai2)]","\n")
  cat("=======================================================================","\n")
  cat("Parameters correlation matrix","\n") 
  print(round(ddd%*%varmat%*%ddd,5))
  llunif<- -n*log(m); csisb<-(m-aver)/(m-1);
  llsb<-loglikcub00(m,freq,1,csisb)
  cat("=======================================================================","\n")
  cat("Log-lik(pai1^,pai2^,csi^)    =",round(loglik,digits=8),"\n")
  cat("Mean Log-likelihood          =",round(loglik/n,digits=8),"\n")
  cat("-----------------------------------------------------------------------","\n")
  cat("Log-lik(UNIFORM)         =",round(llunif,digits=8),"\n")
  cat("Log-lik(Shifted-BINOMIAL)=",round(llsb,digits=8),"\n")
  cat("-----------------------------------------------------------------------","\n")
  theorpr<-probcubshe1(m,pai1,pai2,csi,shelter)
  dissshe<-dissim(theorpr,freq/n)
  llunif<- -n*log(m)
  nonzero<-which(freq!=0)
  logsat<- -n*log(n)+sum((freq[nonzero])*log(freq[nonzero]))
  X2<-sum(((n*theorpr-freq)^2)/(n*theorpr))
  LL2<-1/(1+mean((freq/(n*theorpr)-1)^2))
  II2<-(loglik-llunif)/(logsat-llunif)
  FF2<-1-dissshe
  AICCUBshe<- -2*loglik+2*(3)
  BICCUBshe<- -2*loglik+log(n)*(3)
  ####################################################################
  cat("Pearson Fitting measure    ==>  X^2 =",X2,"(p-val.=",1-pchisq(X2,m-3),")","\n")
  cat("Normed Dissimilarity index ==>  Diss=",round(dissshe,digits=5),"\n")
  cat("F^2 fitting measure        ==>  F^2 =",round(FF2,digits=5),"\n")
  cat("Lik-based fitting measure  ==>  L^2 =",LL2,"\n")
  cat("Relative Log-lik index     ==>  I   =",round(II2,digits=5),"\n")
  cat("-----------------------------------------------------------------------","\n")
  cat("AIC-CUB+she          =",round(AICCUBshe,digits=8),"\n")
  cat("BIC-CUB+she          =",round(BICCUBshe,digits=8),"\n")
  cat("ICOMP-CUB+she        =",round(ICOMP,digits=8),"\n")
  cat("=======================================================================","\n")
  ################################################################
  #        Assignments as global variables
  ################################################################
  #   assign('pai1',pai1,pos=1)
  #   assign('pai2',pai2,pos=1)
  #   assign('csi',csi,pos=1)
  #   assign('varmat',varmat,pos=1)
  #   assign('loglik',loglik,pos=1)
  #   assign('delta',delta,pos=1)
  #   assign('paistar',paistar,pos=1)
  #   assign('nniter',nniter,pos=1)
  ################################################################
  #           /* CUB-she ===>>> Only 1 plot with Diss=... */
  ################################################################
  if (makeplot==TRUE){
    plot(cbind(1:m,1:m),cbind(theorpr,(freq/n)),
         main=paste("CUB model with shelter effect","      (Diss =",round(dissshe,digits=4),")"),
         xlim=c(1,m),ylim=c(0,1.1*max(theorpr,(freq/n))),las=1,
         xlab="Ordinal values of r=1,2,...,m",
         ylab=expression(paste("Observed relative frequencies (dots) and fitted probabilities (circles)")));
    points(1:m,theorpr,pch=21,cex=1.5);
    points(1:m,freq/n,pch=16,cex=1.2);
    ### points(shelter,theorpr[shelter]-delta,pch=8);
    abline(h=0);
  }
  durata<-proc.time()-tt0;durata<-durata[1];
  cat("=======================================================================","\n")
  cat("Elapsed time=",durata,"seconds","=====>>>",date(),"\n")
  cat("=======================================================================","\n","\n")
  results<-list('estimates'=round(stime,digits=5),'loglik'=loglik,'niter'=nniter,'varmat'=varmat,'BIC'=round(BICCUBshe,digits=8))
}
