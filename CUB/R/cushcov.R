# @title CUSH model with covariates
# @description Estimate and validate a CUSH model for ordinal responses, with covariates
#  to explain the shelter effect.
# @usage cushcov(m, ordinal, X, shelter, makeplot)
# @param m Number of ordinal categories
# @param ordinal Vector of ordinal responses
# @param X Matrix of selected covariates for explaining the shelter effect
# @param shelter Category corresponding to the shelter choice
# @param makeplot Logical: if TRUE and if only one dichotomous covariate is included in the model, 
# with levels (0,1), the function returns a graphical plot comparing the distributions of the
#  responses conditioned to the value of the covariate
# @aliases cushcov
# @return An object of the class "CUSH"
# @import stats graphics
#' @keywords internal


cushcov <-
function(m,ordinal,X,shelter,makeplot){ 
  tt0<-proc.time()
  freq<-tabulate(ordinal,nbins=m); n<-length(ordinal);
  fc<-freq[shelter]/n
  delta<-max(0.01,(m*fc-1)/(m-1))              ### sufficient unbiased estimator for a CUSH model
  ncovar<-NCOL(X)
  omzero<-log(delta/(1-delta))         ### initial estimate of omega_0
  omegainit<-c(omzero,rep(0.1,ncovar)) ### initial estimate of omega vector
  ### maximize w.r.t. omega 
  XX<-cbind(1,X)
  esternocush<-cbind(ordinal,XX)
  paravec<-omegainit 
  shelter<-shelter
  optimomega<-optim(paravec,effecush,esternocush,shelter=shelter,m=m,gr=NULL,hessian=TRUE)
  #################################################################
  # Computation of estimates and log-likelihood
  #################################################################
  omegaest<-optimomega$par               #omega estimates
  loglik<-loglikcushcov(m,ordinal,X,omegaest,shelter) #loglik at the maximum
  HHH<-optimomega$hessian
  nparam<-length(omegaest)
  
  if (det(HHH)<=0){
    warning("Variance-Covariance matrix is not positive definite")
    varmat<-ddd<-cormat<-matrix(NA,nrow=nparam,ncol=nparam)
    trvarmat<-ICOMP<-NA
    errst<-wald<-pval<-rep(NA,nparam)  
  } else {
    varmat<-solve(HHH)
    errst<-sqrt(diag(varmat))       ### vector
    ddd<-diag(sqrt(1/diag(varmat))) ### matrix
    wald<-omegaest/errst
    trvarmat<-sum(diag(varmat))
    ICOMP<- -2*loglik + nparam*log(trvarmat/nparam) - log(det(varmat)) ## added
  }
  
  
  AICCUSH<- -2*loglik+2*nparam
  BICCUSH<- -2*loglik+nparam*log(n)
  
  ######## optimomega$message; cat("\n")
  cat("=======================================================================","\n")
  cat(">>>    ML estimation of a CUSH model with covariates  (September 2015)   <<<","\n") 
  cat("=======================================================================","\n")
  cat("n =", n,"    m =",m ,"   shelter =",round(shelter,digits=1), 
      "   Number of covariates for omega =", ncovar,"\n")
  cat("=======================================================================","\n")
  cat("parameters  ML-estimates  stand.errors    Wald-test      p-value ","\n")
  cat("=======================================================================","\n")
  nomi<-c(paste("omega",0:(nparam-1),sep="_"))
  stime<-round(omegaest,5); errstd<-round(errst,5);  wald<-round(wald,5);
  pval<-round(2*(1-pnorm(abs(wald))),20)
  for(i in 1:length(nomi)){
    cat(nomi[i],"     ",stime[i],"      ",errstd[i],"       ",wald[i],"      ",pval[i],"\n")
  }
  cat("=======================================================================","\n")
  cat("Parameters correlation matrix","\n") 
  print(round(ddd%*%varmat%*%ddd,5))
  cat("=======================================================================","\n")
  cat("Log-lik(omega^) =",round(loglik,digits=8),"\n")
  cat("AIC-CUSH-covar  =",round(AICCUSH,digits=8),"\n")
  cat("BIC-CUSH-covar  =",round(BICCUSH,digits=8),"\n")
  cat("ICOMP-CUSH-covar=",round(ICOMP,digits=8),"\n")
  #########################################
  #   assign('varmat',varmat,pos=1)
  #   assign('omega',omegaest,pos=1)
  #   assign('loglik',loglik,pos=1)
  #########################################
  
  
  if(ncovar==1 & length(unique(X))==2) {
    #code dicopai
    vett<-as.matrix(c(0,1))
    delta0<-logis(vett[1],omegaest)
    delta1<-logis(vett[2],omegaest)
    prob0<-probcush(m,delta0,shelter)
    prob1<-probcush(m,delta1,shelter)
    
    maxpr<-max(prob0,prob1)
    if(makeplot==TRUE){
      plot(1:m,prob0,ylim=c(0.0,1.1*maxpr),cex.main=0.8,las=1,
           main="CUSH distributions, given delta-covariate=0, 1",cex=2,
           xlab="",ylab="Prob(R|D=0)  and  Prob(R|D=1)",pch=1,lty=1,type="b");
      lines(1:m,prob1,cex=2,pch=19,lty=2,type="b");
      abline(h=0);
    }
    ### Expected moments given D=0,1 
    exp0<- (1-delta0)*(m+1)/2 + shelter*delta0 ;    exp1<-(1-delta1)*(m+1)/2 + shelter*delta1;
    cushmode0<-which.max(prob0);              cushmode1<-which.max(prob1);
    ### Sample averages and modal value, given D=0,1
    ord0<-ordinal[X==0]; ord1<-ordinal[X==1];
    n0<-length(ord0);    n1<-length(ord1);
    aver0<-mean(ord0);   aver1<-mean(ord1);
    obsmode0<-which.max(tabulate(ordinal[X==0]))    
    obsmode1<-which.max(tabulate(ordinal[X==1]))
    cat("Samples and populations measures, given dichotomous covariate (D=0) and (D=1)","\n")
    cat("-----------------------------------------------------------------------","\n")
    cat("(D = 0)","   n0 = ", n0,
        "        delta_0=",round(delta0,digits=3),"\n")
    cat("............................","\n")
    cat("Sample average  =",round(aver0,digits=8),"   Sample mode =",round(obsmode0,digits=1),"\n")
    cat("CUSH expectation =",round(exp0,digits=8), "    CUSH mode   =",round(cushmode0,digits=1),"\n")
    cat("-----------------------------------------------------------------------","\n")
    cat("(D = 1)","   n1 = ", n1,
        "        delta_1=",round(delta1,digits=3),"\n")
    cat("............................","\n")
    cat("Sample average  =",round(aver1,digits=8)," Sample mode =",round(obsmode1,digits=1),"\n")
    cat("CUSH expectation =",round(exp1,digits=8), " CUSH mode    =",round(cushmode1,digits=1),"\n")
    cat("-----------------------------------------------------------------------","\n")
  }
  
  durata<-proc.time()-tt0;durata<-durata[1];
  cat("=======================================================================","\n")
  cat("Convergence code  =",optimomega$convergence,"\n")
  cat("=======================================================================","\n") 
  cat("Elapsed time      =",durata,"seconds","=====>>>",date(),"\n")
  results<-list('estimates'=stime, 'loglik'=loglik,  'varmat'=varmat,'BIC'=round(BICCUSH,digits=8))
  
}
