# @title Main function for CUBE models with covariates only for feeling
# @description Estimate and validate a CUBE model for ordinal data, with covariates only for explaining the
# feeling component.
# @aliases cubecsi
# @usage cubecsi(m, ordinal, W, starting, maxiter, toler)
# @param m Number of ordinal categories
# @param ordinal Vector of ordinal responses
# @param W Matrix of selected covariates for explaining the feeling component
# @param starting Vector of initial parameters estimates to start the optimization algorithm, with length equal to
#  NCOL(W) + 3 to account for an intercept term for the feeling component (first entry)
# @param maxiter Maximum number of iterations allowed for running the optimization algorithm 
# @param toler Fixed error tolerance for final estimates 
# @return An object of the class "CUBE". For cubecsi, $niter will return a NULL value since the optimization procedure
#  is not iterative but based on "optim" (method = "L-BFGS-B", option hessian=TRUE). \cr $varmat will return the inverse 
#  of the numerically computed hessian when it is positive definite, otherwise the procedure will return a matrix of NA
#   entries.
# @import stats
# @seealso \code{\link{loglikcubecsi}},  \code{\link{inibestcubecsi}},  \code{\link{CUBE}} 
# @examples 
# ### Applying \donttest option since the proposed examples require a long run time for check 
# \donttest{
# data(relgoods)
# m=10
# ordinal=relgoods[,37]
# age=2014-relgoods[,4]
# lage=log(age)-mean(log(age))
# nona=na.omit(cbind(ordinal,lage))
# ordinal=nona[,1]
# W=nona[,2]
# starting=rep(0.1,4)     
# fit=cubecsi(m, ordinal, W, starting, maxiter=100, toler=1e-3)
# param=fit$estimates
# pai=param[1]                        ## ML estimates for the uncertainty parameter
# gama=param[2:(length(param)-1)]     ## ML estimates for the coefficients of the feeling covariates
# phi=param[length(param)]            ## ML estimates for the overdispersion parameter
# loglik=fit$loglik
# varmat=fit$varmat
# BIC=fit$BIC
# ##########################################################
# data(univer)
# m=7 
# ordinal=univer[,8]
# gender=univer[,4]
# initial=inibestcube(m,ordinal)
# starting=inibestcubecsi(m,ordinal,W=gender,initial,maxiter=500,toler=1e-6)
# fitcsi=cubecsi(m, ordinal, W=gender, starting, maxiter=100, toler=1e-3)
# param=fitcsi$estimates
# pai=param[1]                       ## ML estimates for the uncertainty parameter
# gama=param[2:(length(param)-1)]    ## ML estimates for the coefficients of the feeling covariates
# phi=param[length(param)]           ## ML estimates for the overdispersion parameter
# loglik=fitcsi$loglik
# varmat=fitcsi$varmat
# BIC=fitcsi$BIC
# }
#' @keywords internal #models




cubecsi <-
function(m,ordinal,W,starting,maxiter,toler){
  tt0<-proc.time()
  n<-length(ordinal)
  q<-length(starting)-3
  pai<-starting[1]; gama<-starting[2:(q+2)]; phi<-starting[q+3];
  #(0)# log-lik
  loglikzero<-loglikcubecsi(m,ordinal,W,pai,gama,phi)
  #################################################################
  param<-c(pai,gama,phi)
  ##################################
  ### maximize w.r.t. gama and phi #########
  optimparam<-optim(param,effecubecsi,ordinal=ordinal,W=W,m=m,method="L-BFGS-B",lower=c(0.01,rep(-Inf,q+1),0.0001), upper=c(0.99,rep(Inf,q+1),0.5),gr=NULL,hessian=TRUE)
  #################################################################
  # 7.# Computation of updated estimates and log-likelihood
  #################################################################
  paramest<-optimparam$par
  pai<-paramest[1]
  gama<-paramest[2:(q+2)]    #updated gama estimates
  phi<-paramest[q+3]         #updated phi estimates
  hessian<-optimparam$hessian
  ### updated log-likelihood
  loglik<-loglikcubecsi(m,ordinal,W,pai,gama,phi)
  vettestim<-c(pai,gama,phi)
  nparam<-length(vettestim)
  ####################################################################
  AICCUBEcsi<- -2*loglik+2*nparam
  BICCUBEcsi<- -2*loglik+log(n)*nparam
  ###########################################################################
  # Compute asymptotic standard errors of ML estimates via (numerical)Hessian
  ###########################################################################
  
  if (det(hessian)<=0){
    warning("Variance-Covariance matrix is not positive definite")
    varmat<-ddd<-cormat<-matrix(NA,nrow=nparam,ncol=nparam)
    errstd<-wald<-pval<-rep(NA,nparam)
    ICOMP<-trvarmat<-NA
  } else {
    varmat<-solve(hessian)
    errstd<-sqrt(diag(varmat))
    ddd<-diag(sqrt(1/diag(varmat)))
    wald<-vettestim/errstd
    pval<-c(round(2*(1-pnorm(abs(wald))),20))
    cormat<-(ddd%*%varmat)%*%ddd
    trvarmat<-sum(diag(varmat))
    ICOMP<- -2*loglik + nparam*log(trvarmat/nparam) - log(det(varmat))
    
    errstd<-round(errstd,5)  
    wald<-round(wald,5)
    pval<-round(pval,5)
  }
  
  nomi<-c("pai    ",paste("gamma",0:(length(gama)-1),sep="_"),"phi    ")
  stime<-round(vettestim,5)
  ####################################################################
  ### Print CUBEcsi results of ML estimation  
  ####################################################################
  cat("\n")
  cat("=======================================================================","\n")
  cat("==> CUBEcsi Program <<<=== ML-estimates via optim post E-M algorithm   ","\n")
  cat("=======================================================================","\n")
  cat("                Covariates for feeling ==> q =", q,"\n")
  cat("=======================================================================","\n")
  cat("  *** m=", m,"       *** Sample size: n=", n,"                         ", "\n")
  cat("=======================================================================","\n")
  cat("parameters  ML-estimates  stand.errors    Wald-test      p-value ","\n")
  cat("=======================================================================","\n")
  for(i in 1:length(nomi)){
    cat(nomi[i],"     ",stime[i],"      ",errstd[i],"       ",wald[i],"      ",pval[i],"\n")
  }
  ####################################################################
  cat("=======================================================================","\n")
  cat("                         Parameters correlation matrix","\n") 
  rownames(cormat)<-nomi; colnames(cormat)<-nomi; 
  print(round(cormat,3))
  ##############################################################################
  cat("=======================================================================","\n")
  cat("Log-lik(pai^,gama^,phi^) =",round(loglik,digits=8),"\n")
  cat("Mean Log-likelihood      =",round(loglik/n,digits=8),"\n")
  cat("-----------------------------------------------------------------------","\n")
  cat("AIC-CUBE-csi       =",round(AICCUBEcsi,digits=8),"\n")
  cat("BIC-CUBECSI        =",round(BICCUBEcsi,digits=8),"\n")
  cat("ICOMP-CUBECSI      =",round(ICOMP,digits=8),"\n")
  cat("=======================================================================","\n")  
  ################################################################
  #        Assignments as global variables
  ################################################################
  #   assign('pai',pai,pos=1)
  #   assign('gama',gama,pos=1)
  #   assign('phi',phi,pos=1)
  #   assign('varmat',varmat,pos=1)
  #   assign('loglik',loglik,pos=1)
  durata<-proc.time()-tt0;durata<-durata[1];
  cat("=======================================================================","\n")
  cat("Convergence code =",optimparam$convergence,"\n")
  cat("=======================================================================","\n")  
  cat("Elapsed time     =",durata,"seconds","=====>>>",date(),"\n")
  results<-list('estimates'=stime, 'loglik'=loglik, 'varmat'=varmat,'BIC'= round(BICCUBEcsi,digits=8))
}
