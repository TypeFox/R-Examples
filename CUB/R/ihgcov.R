# @title Main function for IHG models with covariates
# @description Estimate and validate an IHG model for given ordinal responses, with covariates to 
# explain the preference parameter. 
# @aliases ihgcov
# @usage ihgcov(m, ordinal, U, makeplot)
# @param m Number of ordinal categories
# @param ordinal Vector of ordinal responses
# @param U Matrix of selected covariates for the preference parameter 
# @param makeplot Logical: if TRUE and if only one dichotomous covariate is included in the model, with levels (0,1), 
# the function returns a graphical plot comparing the distributions of the responses conditioned to the value of the 
# covariate
# @return An object of the class "IHG"
# @details The optimization procedure is run via "optim", option method="Brent" for constrained optimization 
# (lower bound = 0, upper bound=1).
# @import stats graphics
# @return An object of the class "IHG"
#' @keywords internal 


ihgcov <-
function(m,ordinal,U,makeplot){ 
  tt0<-proc.time()
  freq<-tabulate(ordinal,nbins=m); n<-length(ordinal);
  theta<-iniihg(m,freq)
  ncovar<-NCOL(U)
  nuzero<-log(theta/(1-theta))         ### initial estimate of nu_0
  nuinit<-c(nuzero,rep(0.1,ncovar)) ### initial estimate of nu vector
  ### maximize w.r.t. nu
  nu<-nuinit
  optimnu<-optim(nu,effeihgcov,ordinal=ordinal,U=U,m=m,hessian=TRUE)
  #################################################################
  # Computation of estimates and log-likelihood
  #################################################################
  nuest<-optimnu$par                #nu estimates
  loglik<-loglikihgcov(m,ordinal,U,nuest)
  HHH<-optimnu$hessian
  nparam<-length(nuest)
  
  if (det(HHH)<=0){
    warning("Variance-covariance matrix not-positive definite")
    varmat<-ddd<-matrix(NA,nrow=nparam,ncol=nparam)
    errst<-wald<-rep(NA,nparam)
    trvarmat<-ICOMP<-NA    
  } else {
    varmat<-solve(HHH)
    errst<-sqrt(diag(varmat))       ### vector
    ddd<-diag(sqrt(1/diag(varmat))) ### matrix
    wald<-nuest/errst
    trvarmat<-sum(diag(varmat))
    ICOMP<- -2*loglik + nparam*log(trvarmat/nparam) - log(det(varmat)) 
  }
  
  AICIHGCOV<- -2*loglik+2*nparam
  BICIHGCOV<- -2*loglik+nparam*log(n)
  cat("=======================================================================","\n")
  cat(">>>    ML estimation of an IHG model with covariates  (September 2015)    <<<","\n") 
  cat("=======================================================================","\n")
  cat("n =", n,"    m =",m ,"  
      Number of covariates for nu =", ncovar,"\n")
  cat("=======================================================================","\n")
  cat("parameters  ML-estimates  stand.errors    Wald-test      p-value ","\n")
  cat("=======================================================================","\n")
  nomi<-c(paste("nu",0:(nparam-1),sep="_"))
  stime<-round(nuest,5); errstd<-round(errst,5);  wald<-round(wald,5);
  pval<-round(2*(1-pnorm(abs(wald))),20)
  for(i in 1:length(nomi)){
    cat(nomi[i],"     ",stime[i],"      ",errstd[i],"       ",wald[i],"      ",pval[i],"\n")
  }
  cat("=======================================================================","\n")
  cat("Parameters correlation matrix","\n") 
  print(round(ddd%*%varmat%*%ddd,5))
  cat("=======================================================================","\n")
  cat("Log-lik(nu^)   =",round(loglik,digits=8),"\n")
  cat("AIC-IHG-covar  =",round(AICIHGCOV,digits=8),"\n")
  cat("BIC-IHG-covar  =",round(BICIHGCOV,digits=8),"\n")
  cat("ICOMP-IHG-covar=",round(ICOMP,digits=8),"\n")
  cat("=======================================================================","\n")
  #########################################
  #   assign('varmat',varmat,pos=1)
  #   assign('nu',nuest,pos=1)
  #   assign('loglik',loglik,pos=1)
  #########################################
  
  if(ncovar==1 & length(unique(U))==2) {
    #code dicopai
    vett<-as.matrix(c(0,1))
    theta0<-logis(vett[1],nuest)
    theta1<-logis(vett[2],nuest)
    prob0<-probihg(m,theta0)
    prob1<-probihg(m,theta1)
    maxpr<-max(prob0,prob1)
    if(makeplot==TRUE){
      plot(1:m,prob0,ylim=c(0.0,1.1*maxpr),cex.main=0.8,las=1,
           main="IHG distributions, given theta-covariate=0, 1",cex=2,
           xlab="",ylab="Prob(R|D=0)  and  Prob(R|D=1)",pch=1,lty=1,type="b");
      lines(1:m,prob1,cex=2,pch=19,lty=2,type="b");
      abline(h=0);
    }
    ### Expected moments given D=0,1 
    exp0<-(m-theta0)/(1+(theta0)*(m-2)) ;     exp1<-(m-theta1)/(1+(theta1)*(m-2));
    ihgmode0<-which.max(prob0);              ihgmode1<-which.max(prob1);
    ### Sample averages and modal value, given D<-0,1
    ord0<-ordinal[U==0]; ord1<-ordinal[U==1];
    n0<-length(ord0);    n1<-length(ord1);
    aver0<-mean(ord0);   aver1<-mean(ord1);
    obsmode0<-which.max(tabulate(ordinal[U==0]))    
    obsmode1<-which.max(tabulate(ordinal[U==1]))
    cat("Samples and populations measures, given dichotomous covariate (D=0) and (D=1)","\n")
    cat("-----------------------------------------------------------------------","\n")
    cat("(D = 0)","   n0 = ", n0,
        "        theta_0=",round(theta0,digits=3),"\n")
    cat("............................","\n")
    cat("Sample average  =",round(aver0,digits=8),"   Sample mode =",round(obsmode0,digits=1),"\n")
    cat("IHG expectation =",round(exp0,digits=8), "    IHG mode   =",round(ihgmode0,digits=1),"\n")
    cat("-----------------------------------------------------------------------","\n")
    cat("(D = 1)","   n1 = ", n1,
        "        theta_1=",round(theta1,digits=3),"\n")
    cat("............................","\n")
    cat("Sample average  =",round(aver1,digits=8)," Sample mode =",round(obsmode1,digits=1),"\n")
    cat("IHG expectation =",round(exp1,digits=8), " IHG mode    =",round(ihgmode1,digits=1),"\n")
    cat("-----------------------------------------------------------------------","\n")
  }
  
  
  
  durata<-proc.time()-tt0;durata<-durata[1];
  cat("=======================================================================","\n")
  cat("Convergence code =",optimnu$convergence,"\n")
  cat("=======================================================================","\n")
  cat("Elapsed time =",durata,"seconds","=====>>>",date(),"\n")
  results<-list('estimates'=stime, 'loglik'=loglik, 'varmat'=varmat,'BIC'=round(BICIHGCOV,digits=8))
  
}
