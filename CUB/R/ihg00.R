# @title Main function for IHG models without covariates
# @description Estimate and validate an IHG model without covariates for given ordinal responses. 
# @aliases ihg00
# @usage ihg00(m, ordinal, makeplot)
# @param m Number of ordinal categories
# @param ordinal Vector of ordinal responses
# @param makeplot Logical: if TRUE, the function returns a graphical plot comparing fitted probabilities
#  and observed relative frequencies
# @return An object of the class "IHG"
# @details The optimization procedure is run via "optim", option method="Brent" for constrained optimization 
# (lower bound = 0, upper bound=1).
# @import stats
# @return An object of the class "IHG"
#' @keywords internal 



ihg00 <-
function(m,ordinal,makeplot){
  tt0<-proc.time()
  freq<-tabulate(ordinal,nbins=m)
  n<-sum(freq)
  aver<-mean(ordinal)
  theta<-iniihg(m,freq) ### initial value (moment estimator of theta)
  estoptim<-optim(theta,effeihg,m=m,freq=freq,method="Brent",lower=0, upper=1,hessian=TRUE)
  theta<-estoptim$par
  errstdtheta<-1/sqrt(estoptim$hessian)
  varmat<-errstdtheta^2
  loglik<-loglikihg(m,freq,theta)
  wald2<-(theta-1/m)/errstdtheta
  pvaltheta2<-round(2*(1-pnorm(abs(wald2))),20)
  ##wald=theta/errstdtheta
  # prima: pvaltheta1=round(2*(1-pnorm(abs(wald))),20)
  # prima: pvaltheta1=round(2*(1-pnorm(abs(wald))),20)
  AICIHG<- -2*loglik+2
  BICIHG<- -2*loglik+log(n)
  #####################################################
  cat("\n")
  cat("=======================================================================","\n")
  cat(">>>    ML estimation of an IHG model without covariates  (September 2015) <<","\n") 
  cat("=======================================================================","\n")
  cat("*** m=", m,"                *** Sample size: n=", n,"\n")
  cat("=======================================================================","\n")
  cat("parameter  ML-estimates  stand.errors   Wald-test  Test vs H0:theta=1/m ","\n")
  cat("=======================================================================","\n")
  cat("theta      ", theta,"   ",errstdtheta,"   ",wald2,"      ",pvaltheta2,"\n")
  cat("=======================================================================","\n")
  csisb<-(m-aver)/(m-1); llsb<-loglikcub00(m,freq,1,csisb);
  llunif<- -n*log(m)
  nonzero<-which(freq!=0)
  logsat<- -n*log(n)+sum((freq[nonzero])*log(freq[nonzero]))
  cat("Log-lik(theta^)    =",round(loglik,digits=8),"\n")
  cat("Mean Log-likelihood=",round(loglik/n,digits=8),"\n")
  cat("-----------------------------------------------------------------------","\n")
  cat("Log-lik(UNIFORM)         =",round(llunif,digits=8),"\n")
  cat("Log-lik(Shifted-BINOMIAL)=",round(llsb,digits=8),"\n")
  cat("Log-lik(Saturated)       =",round(logsat,digits=8),"\n")
  cat("-----------------------------------------------------------------------","\n")
  theorpr<-probihg(m,theta)
  dissihg<-dissim(theorpr,freq/n)
  X2<-sum(((n*theorpr-freq)^2)/(n*theorpr))
  LL2<-1/(1+mean((freq/(n*theorpr)-1)^2))
  II2<-(loglik-llunif)/(logsat-llunif)
  FF2<-1-dissihg
  cat("Pearson Fitting measure    ==>  X^2 =",X2,"(p-val.=",1-pchisq(X2,m-3),")","\n")
  cat("Normed Dissimilarity index ==>  Diss=",round(dissihg,digits=5),"\n")
  cat("F^2 fitting measure        ==>  F^2 =",round(FF2,digits=5),"\n")
  cat("Lik-based fitting measure  ==>  L^2 =",LL2,"\n")
  cat("Relative Log-lik index     ==>  I   =",round(II2,digits=5),"\n")
  cat("=======================================================================","\n")
  cat("AIC-IHG         =",round(AICIHG,digits=8),"\n")
  cat("BIC-IHG         =",round(BICIHG,digits=8),"\n")
  ################################################################
  #        Assignments as global variables
  ################################################################
  #   assign('theta',round(theta,digits=5),pos=1)
  #   assign('loglik',loglik,pos=1)
  #   assign('varmat',varmat,pos=1)
  ################################################################
  #           /* IHG ===>>> Only 1 plot with Diss=... */
  ################################################################
  if (makeplot==TRUE){
    plot(cbind(1:m,1:m),cbind(theorpr,(freq/n)),
         main=paste("IHG model (without covariates)","      (Diss =",round(dissihg,digits=4),")"),
         xlim=c(1,m),ylim=c(0,1.1*max(theorpr,(freq/n))),las=1,
         xlab="Ordinal values of r=1,2,...,m",
         ylab=expression(paste("Observed relative frequencies (dots) and fitted probabilities (circles)")))
    points(1:m,theorpr,pch=21,cex=1.5)
    points(1:m,freq/n,pch=16,cex=1.2)
    ### points(shelter,theorpr[shelter]-delta,pch=8);
    abline(h=0);
  }
  durata<-proc.time()-tt0;durata<-durata[1];
  cat("=======================================================================","\n")
  cat("Convergence code  =",estoptim$convergence,"\n")
  cat("=======================================================================","\n")
  cat("Elapsed time      =",durata,"seconds","=====>>>",date(),"\n")
  stime=round(theta,digits=5)
  results<-list('estimates'=stime, 'loglik'=loglik,'varmat'=round(varmat,digits=5),'BIC'=round(BICIHG,digits=8))
  
}
