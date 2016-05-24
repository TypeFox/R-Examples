# @title Main function for CUBE models with covariates
# @description Function to estimate and validate a CUBE model with 
# explicative covariates for all the three parameters.
# @aliases cubecov
# @usage cubecov(m, ordinal, Y, W, Z, starting, maxiter, toler)
# @param m Number of ordinal categories
# @param ordinal Vector of ordinal responses
# @param Y Matrix of selected covariates for explaining the uncertainty component
# @param W Matrix of selected covariates for explaining the feeling component
# @param Z Matrix of selected covariates for explaining the overdispersion component
# @param starting Vector of initial parameters estimates to start the optimization algorithm 
# (it has length NCOL(Y) + NCOL(W) + NCOL(Z) + 3 to account for intercept terms 
# for all the three components
# @param maxiter Maximum number of iterations allowed for running the optimization algorithm 
# @param toler Fixed error tolerance for final estimates 
# @return An object of the class "CUBE"
# @import stats
# @references
# Piccolo, D. (2014). Inferential issues on CUBE models with covariates,
#  \emph{Communications in Statistics - Theory and Methods}, \bold{44}, 
#  DOI: 10.1080/03610926.2013.821487
#' @keywords internal


cubecov <-
function(m,ordinal,Y,W,Z,starting,maxiter,toler){
  tt0<-proc.time()
  n<-length(ordinal)
  p<-NCOL(Y)
  q<-NCOL(W)
  v<-NCOL(Z)
  aver<-mean(ordinal)
  YY<-cbind(1,Y)   
  WW<-cbind(1,W)       
  ZZ<-cbind(1,Z)    
  # **********************************
  # *** E-M algorithm for CUBECOV  ***
  # **********************************
  #################################################################
  # 00.# Initial values of parameters ............. Attention !!! #
  #################################################################
  betjj<-starting[1:(p+1)]
  gamajj<-starting[(p+2):(p+q+2)]
  alphajj<-starting[(p+q+3):(p+q+v+3)]
  ################################
  # * * * Iterative Loop    * * *#
  ################################
  nniter<-1
  while(nniter<=maxiter){
    #################################################################
    # 01.# Initial values of vectors i=1,2,..,n
    #################################################################
    paijj<-logis(Y,betjj); csijj<-logis(W,gamajj) 
    phijj<-1/(1/logis(Z,alphajj)-1)  #****
    #################################################################
    # 02.# Computation of beta-binomial distribution i=1,2,..,n
    #################################################################
    betabin<-betabinomial(m,ordinal,csijj,phijj)
    #################################################################
    # 03.# Computation of CUBE probability distribution i=1,2,..,n
    #################################################################
    probi<-paijj*(betabin-1/m)+1/m
    likold<-sum(log(probi))
    #################################################################
    # 4.# Computation of conditional probability i=1,2,..,n
    #################################################################
    taui<-1-(1-paijj)/(m*probi)
    #################################################################
    # 5. Unify parameter vectors
    param<-c(gamajj,alphajj)
    #################################################################
    # 6.# Maximization of Q_1(beta) and Q_2(param)
    #################################################################
    ### maximize w.r.t. bet and gama #########
    esterno1<-cbind(taui,YY) 
    covar<-esterno1[,2:NCOL(esterno1)]
    esterno2<-cbind(taui,ordinal,W,Z) 
    bet<-betjj
    optimbet<-optim(bet,Quno,esterno1=esterno1,gr=NULL) #added gr  
    optimparam<-optim(param,Qdue,esterno2=esterno2,q=q,m=m,gr=NULL)
    #################################################################
    # 7.# Computation of updated estimates and log-likelihood
    #################################################################
    betjj<-optimbet$par                #updated bet estimates
    paramjj<-optimparam$par
    gamajj<-paramjj[1:(q+1)]           #updated gama estimates
    alphajj<-paramjj[(q+2):(q+v+2)]    #updated alpha estimates
    ### updated log-likelihood
    liknew<-loglikcubecov(m,ordinal,Y,W,Z,betjj,gamajj,alphajj)
    #################################################################
    # 8.# Checking improvement of updated log-likelihood
    #################################################################
    # print(c(nniter,betjj,gamajj,loglikjj)); #OPTIONAL PRINTING OF ITERATIONS
    testloglik<-abs(liknew-likold)
    #print(nniter);##added
    #print(round(c(liknew,likold,testloglik),digits=7))#added
    if(testloglik<=toler) break else {likold<-liknew}
    nniter<-nniter+1
  }
  #################################################################
  # 8.# Final ML estimates and maximized log-likelihood
  #################################################################
  bet<-betjj;  gama<-gamajj;  alpha<-alphajj;
  loglik<-liknew
  ###### End of E-M algorithm for CUBE ***********************************************
  paramest<-c(bet,gama,alpha)
  nparam<- length(paramest) ###p+q+v+3;
  AICCUBE<- -2*loglik+2*nparam
  BICCUBE<- -2*loglik+log(n)*nparam
  
  ############################################################
  # Compute asymptotic standard errors of ML CUBE estimates ##
  ############################################################
  
  varmat<-varcovcubecov(m,ordinal,Y,W,Z,bet,gama,alpha)
  #if(det(varmat)<=0) stop("Variance-Covariance matrix NOT positive definite")
  nomi<-c(paste("beta",0:(length(bet)-1),sep="_"),
         paste("gamma",0:(length(gama)-1),sep="_"),
         paste("alpha",0:(length(alpha)-1),sep="_"))
  stime<-round(paramest,5)
  
  if (isTRUE(varmat==matrix(NA,nrow=nparam,ncol=nparam))==TRUE){
    ddd<-cormat<-matrix(NA,nrow=nparam,ncol=nparam)
    trvarmat<-ICOMP<-NA
    errstd<-wald<-pval<-rep(NA,nparam)
  } else {
    ddd<-diag(sqrt(1/diag(varmat)))
    cormat<-(ddd%*%varmat)%*%ddd
    trvarmat<-sum(diag(varmat))
    ICOMP<- -2*loglik + nparam*log(trvarmat/nparam) - log(det(varmat))
    errstd<-round(sqrt(diag(varmat)),5);  wald<-round(stime/errstd,5);
    pval<-round(2*(1-pnorm(abs(wald))),20)
  }
  
  ##################################################
  # Print CUBE-covariates results of ML estimation #
  ##################################################
  cat("\n")
  cat("=======================================================================","\n")
  cat("==============  CUB Program: version 4.0 (September 2015) =============","\n")
  cat("=======================================================================","\n")
  cat(">> C U B E model with covariates <<   ML-estimates via E-M algorithm   ","\n")
  cat("=======================================================================","\n")
  cat("Covar.pai: p=", p,";  Covar.csi: q=", q,";  Covar.phi: v=",v,  "\n")
  cat("=======================================================================","\n")
  cat("*** m=", m,"  *** Sample size: n=", n,"   *** Iterations=",nniter,"Maxiter=",maxiter,"\n")
  cat("=======================================================================","\n")
  cat("parameters  ML-estimates  stand.errors    Wald-test      p-value ","\n")
  cat("=======================================================================","\n")
  for(i in 1:length(nomi)){
    cat(nomi[i],"     ",stime[i],"      ",errstd[i],"       ",wald[i],"      ",pval[i],"\n")
  }
  ####################################################################
  cat("=======================================================================","\n")
  cat("                         Parameters correlation matrix","\n") 
  rownames(cormat)<-nomi;colnames(cormat)<-nomi; 
  print(round(cormat,5))
  ##############################################################################
  cat("=======================================================================","\n")
  cat("Log-lik(beta^,gamma^,alpha^)=",round(loglik,digits=8),"\n")
  cat("Mean Log-likelihood         =",round(loglik/n,digits=8),"\n")
  cat("-----------------------------------------------------------------------","\n")
  cat("AIC-CUBE           =",round(AICCUBE,digits=8),"\n")
  cat("BIC-CUBE           =",round(BICCUBE,digits=8),"\n")
  cat("ICOMP-CUBE         =",round(ICOMP,digits=8),"\n")
  cat("=======================================================================","\n")
  #####################################################################
  # Assignments as global variables: assign('name',value,pos=1)
  #####################################################################
  #   assign('nniter',nniter,pos=1)
  #   assign('bet',bet,pos=1)
  #   assign('gama',gama,pos=1)
  #   assign('alpha',alpha,pos=1)
  #   assign('loglik',loglik,pos=1)
  #   assign('varmat',varmat,pos=1)
  ####################################
  durata<-proc.time()-tt0;durata<-durata[1];
  cat("=======================================================================","\n")
  cat("Convergence code =",optimparam$convergence,"\n")
  cat("=======================================================================","\n")
  cat("Elapsed time     =",durata,"seconds","=====>>>",date(),"\n")
  
  results<-list('estimates'=stime, 'loglik'=loglik, 'niter'= nniter, 'varmat'=varmat,'BIC'=round(BICCUBE,digits=8))
  
}
