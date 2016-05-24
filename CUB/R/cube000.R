# @title Main function for CUBE models without covariates
# @description Estimate and validate a CUBE model without covariates.
# @aliases cube000
# @usage cube000(m, ordinal, starting, maxiter, toler, makeplot, expinform)
# @param m Number of ordinal categories
# @param ordinal Vector of ordinal responses
# @param starting Vector of initial estimates to start the optimization algorithm, 
# whose length equals the number of parameters of the model
# @param maxiter Maximum number of iterations allowed for running the optimization algorithm 
# @param toler Fixed error tolerance for final estimates 
# @param makeplot Logical: if TRUE, the function returns a graphical plot comparing fitted probabilities 
# and observed relative frequencies
# @param expinform Logical: if TRUE, the function returns the expected variance-covariance matrix
# @return An object of the class "CUBE"
# @import stats
# @references
#  Iannario, M. (2014). Modelling Uncertainty and Overdispersion in Ordinal Data,
#   \emph{Communications in Statistics - Theory and Methods}, \bold{43}, 771--786 \cr
# Iannario, M. (2015). Detecting latent components in ordinal data with overdispersion by means 
# of a mixture distribution, \emph{Quality & Quantity}, \bold{49}, 977--987
# @examples 
# ### Applying donttest option since the proposed examples require long time run for check 
# \donttest{
# data(relgoods)
# m=10
# ordinal=na.omit(relgoods[,37])
# starting = rep(0.1, 3)                              
# fitcube=cube000(m, ordinal, starting, maxiter=500, toler=1e-6, makeplot=TRUE, expinform=FALSE)
# param=fitcube$estimates
# pai=param[1]           # ML estimate for the uncertainty parameter
# csi=param[2]           # ML estimate for the feeling parameter
# phi=param[3]           # ML estimate for the overdispersion parameter
# maxlik=fitcube$loglik 
# niter=fitcube$niter
# BIC=fitcube$BIC
# ###################
# data(univer)
# m=7
# ordinal=univer[,8]
# starting=inibestcube(m,ordinal)    
# model=cube000(m,ordinal,starting,maxiter=200,toler=1e-4,makeplot=TRUE,expinform=TRUE)
# param=model$estimates   # Final ML estimates (pai,csi,phi)
# maxlik=model$loglik
# model$varmat
# model$niter
# model$BIC
# }
#' @keywords internal


cube000 <-
function(m,ordinal,starting,maxiter,
                  toler,makeplot,expinform){ #default for expinform = FALSE
  tt0<-proc.time()
  freq<-tabulate(ordinal,nbins=m); n<-sum(freq); 
  aver<-mean(ordinal); varcamp<-mean(ordinal^2)-aver^2;
  ########################################################
  #(00)# initial estimates, not efficient:
  starting<-inibestcube(m,ordinal)
  pai<-starting[1]; csi<-starting[2]; phi<-starting[3];
  #(0)# log-lik
  loglik<-loglikcube(m,freq,pai,csi,phi)
  
  # ********************************************************************
  # *************   E-M algorithm for CUBE     *************************
  # ********************************************************************
  nniter<-1
  while(nniter<=maxiter){
    likold<-loglik
    
    #(1)# betar
    bb<-betar(m,csi,phi)
    aa<-(1-pai)/(m*pai*bb)
    
    #(2)# taunor
    tauno<-1/(1+aa)
    
    #(3)# pai(k+1)
    pai<-sum(freq*tauno)/n # updated pai estimate
    
    paravecjj<-c(csi,phi)
    
    #(4)# Q(k+1)
    dati<-cbind(tauno,freq)
    ################ EFFECUBE is Q(csi,phi) ###########################
    #(5)# (csi(k+1),phi(k+1))
    ################################## maximize w.r.t. paravec  ########
    paravec<-paravecjj
    optimestim<-optim(paravec,effecube,dati=dati,m=m,method = "L-BFGS-B",lower=c(0.01,0.01),upper=c(0.99,0.5))    # print(nlmaxg)
    ################################################################         
    
    #(6)# theta(k+1)
    paravecjj<-optimestim$par   # updated paravec estimates
    csi<-paravecjj[1];   phi<-paravecjj[2];
    ##########################################
    if(pai<0.001){pai<-0.001; nniter<-maxiter-1}
    #if(csi<0.001){csi<-0.001; nniter<-maxiter-1}
    #if(phi<0.001){phi<-0.001; nniter<-maxiter-1}
    if(pai>0.999){pai<-0.99}         ### to avoid division by 0 !!!
    #if(csi>0.999){csi<-0.99}         ### to avoid division by 0 !!!
    ###################################### print(c(nniter,pai,csi,phi));
    
    #(7)# elle(theta(k+1))
    liknew<-loglikcube(m,freq,pai,csi,phi)
    
    #(8)# test
    testll<-abs(liknew-likold)           # OPTIONAL printing: print(testll); 
    # OPTIONAL printing: print(cbind(nniter,testll,pai,csi,phi));
    if(testll<=toler) break else {loglik<-liknew} # OPTIONAL printing: print(loglik);
    nniter<-nniter+1
  }
  loglik<-liknew
  ###### End of E-M algorithm for CUBE ***********************************************
  AICCUBE<- -2*loglik+2*(3)
  BICCUBE<- -2*loglik+log(n)*(3)
  # ********************************************************
  # Compute ML var-cov matrix and print result for CUBE
  # ********************************************************
  
  if(expinform==TRUE){
    varmat<-varcovcubeexp(m,pai,csi,phi,n)
  }
  else{
    varmat<-varcovcubeobs(m,pai,csi,phi,freq)
  }
  
  nomi<-rbind("pai  ","csi  ","phi  ")
  stime<-c(pai,csi,phi)
  nparam<-length(stime)
  if (isTRUE(varmat==matrix(NA,nrow=nparam,ncol=nparam))==TRUE){
    ddd<-matrix(NA,nrow=nparam,ncol=nparam)
    trvarmat<-ICOMP<-NA
    errstd<-wald<-pval<-rep(NA,nparam)
  } else {
    trvarmat<-sum(diag(varmat))
    ICOMP<- -2*loglik + nparam*log(trvarmat/nparam) - log(det(varmat)) 
    errstd<-sqrt(diag(varmat))
    print(errstd)
    wald<-stime/errstd
    pval<-round(2*(1-pnorm(abs(wald))),20)
    ddd<-diag(sqrt(1/diag(varmat)))
  }
  ####################################################################
  # Print CUBE results of ML estimation  
  ####################################################################
  cat("\n")
  cat("=======================================================================","\n")
  cat("============== CUB Program: version 4.0 (September 2015) ==============","\n")
  cat("=======================================================================","\n")
  cat("> CUBE model Inference                                                 ","\n")
  cat("=======================================================================","\n")
  if(expinform==TRUE){
    cat("Maximum Likelihood estimates (E-M algorithm) ..... Expected Information","\n")
  }
  else{
    cat("Maximum Likelihood estimates (E-M algorithm) ..... Observed Information","\n")
  }
  cat("=======================================================================","\n")
  cat("*** m=", m," *** Sample size: n=", n,"  *** Iterations=",nniter,"Maxiter=",maxiter,"\n")
  cat("=======================================================================","\n")
  cat("parameters  ML-estimates    stand.errors    estimates/stand.errors       ","\n")
  cat("=======================================================================","\n")
  for(i in 1:3){
    cat(nomi[i],"        ",round(stime[i],5),"        ",round(errstd[i],5),"        ",round(wald[i],5),"      ","\n")
  }
  cat("=======================================================================","\n")
  cat("Parameters correlation matrix","\n") 
  print(round(ddd%*%varmat%*%ddd,5))
  ### Log-likelihood comparisons ##############################################
  llunif<- -n*log(m) 
  csisb<-(m-aver)/(m-1)
  llsb<-loglikcub00(m,freq,1,csisb)
  nonzero<-which(freq!=0)
  logsat<- -n*log(n)+sum((freq[nonzero])*log(freq[nonzero]))
  devian<-2*(logsat-loglik)
  cat("\n")
  cat("==================================================== Log-likelihoods ===","\n")
  cat("Log-lik(pai^,csi^,phi^)  =",round(loglik,digits=8),"\n")
  cat("Log-lik(Saturated)       =",round(logsat,digits=8),"\n")
  cat("Log-lik(Shifted-Binomial)=",round(llsb,digits=8),"\n")
  cat("Log-lik(Uniform)         =",round(llunif,digits=8),"\n")
  cat("-----------------------------------------------------------------------","\n")
  cat("Mean Log-likelihood=",round(loglik/n,digits=8),"\n")
  cat("Deviance           =",round(devian,digits=8),"\n")
  ### Fitting measures #########################################################
  cat("=========================================== Global fitting measures ===","\n")
  cat("\n")
  theorpr<-probcube(m,pai,csi,phi)
  dissimcube<-dissim(theorpr,freq/n)
  pearson<-((freq-n*theorpr))/sqrt(n*theorpr)
  X2<-sum(pearson^2)
  relares<-(freq/n-theorpr)/theorpr
  FF2<-1-dissimcube
  LL2<-1/(1+mean((freq/(n*theorpr)-1)^2))
  II2<-(loglik-llunif)/(logsat-llunif)
  cat("Pearson Fitting measure    ==>  X^2 =",X2,"(p-val.=",1-pchisq(X2,m-4),")","\n")
  cat("F^2 fitting measure        ==>  F^2 =",round(FF2,digits=5),"\n")
  cat("Normed Dissimilarity index ==>  Diss=",round(dissimcube,digits=5),"\n")
  cat("Lik-based fitting measure  ==>  L^2 =",round(LL2,digits=5),"\n")
  cat("Relative Log-lik index     ==>  I   =",round(II2,digits=5),"\n")
  cat("-----------------------------------------------------------------------","\n")
  cat("AIC-CUBE          =",round(AICCUBE,digits=8),"\n")
  cat("BIC-CUBE          =",round(BICCUBE,digits=8),"\n")
  cat("ICOMP-CUBE        =",round(ICOMP,digits=8),"\n")
  cat("=======================================================================","\n")
  cat("(R=r) Observed CUBE-prob","Pearson","Relative res.","\n")
  stampa<-cbind(1:m,freq/n,theorpr,pearson,relares)
  print(stampa,digits=5)
  #############################################################################
  #   /* CUBE ===> plot with Diss=.......... optioned by makeplot=TRUE/FALSE */
  #############################################################################
  if(makeplot==TRUE){
    stringtitle="CUBE model estimation ";
    plot(cbind(1:m,1:m),cbind(theorpr,(freq/n)),las=1,
         main=paste(stringtitle,  "     (Diss =",round(dissimcube,digits=4),")"),
         xlim=c(1,m),ylim=c(0.0,1.1*max(theorpr,(freq/n))),
         xlab="Ordinal values of R=1,2,...,m",
         ylab="Observed relative frequencies (dots) and fitted probabilities (circles)");
    ###
    points(1:m,theorpr,pch=21,cex=1.5,lwd=2.0,type="b",lty=3); ### ex pch=8,col="red"
    points(1:m,freq/n,pch=16,cex=1.25,lwd=1.5);
    abline(h=0);
  }
  #####################################################################
  # Assignments as global variables: assign('name',value,pos=1)
  #####################################################################
  #   assign('nniter',nniter,pos=1)
  #   assign('pai',pai,pos=1)
  #   assign('csi',csi,pos=1)
  #   assign('phi',phi,pos=1)
  #   assign('varmat',varmat,pos=1)
  #   assign('loglik',loglik,pos=1)
  ####################################
  durata<-proc.time()-tt0;durata<-durata[1];
  cat("=======================================================================","\n")
  cat("Convergence code =",optimestim$convergence,"\n")
  cat("=======================================================================","\n")
  cat("Elapsed time     =",durata,"seconds","=====>>>",date(),"\n")
  results<-list('estimates'= round(stime, digits=5), 'loglik'= loglik, 'niter'= nniter, 'varmat'= varmat,'BIC'=round(BICCUBE,digits=8))
  
}
