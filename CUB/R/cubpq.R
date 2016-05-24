# @title Main function for CUB models with covariates for both the uncertainty and the feeling components
# @description Estimate and validate a CUB model for given ordinal responses, with covariates for explaining both the
#  feeling and the uncertainty components by means of logistic transform.
# @aliases cubpq
# @usage cubpq(m, ordinal, Y, W, maxiter, toler, makeplot)
# @param m Number of ordinal categories
# @param ordinal Vector of ordinal responses
# @param Y Matrix of selected covariates for explaining the uncertainty component
# @param W Matrix of selected covariates for explaining the feeling component
# @param maxiter Maximum number of iterations allowed for running the optimization algorithm 
# @param toler Fixed error tolerance for final estimates 
# @param makeplot Logical: if TRUE, and if only a dichotomous covariate is included in the model to 
# explain both the uncertainty and the feeling components, with levels (0,1),  the function returns a graphical plot
#  comparing the distributions of the responses conditioned to the value of the covariate
# @return An object of the class "CUB"
# @import stats
# @seealso \code{\link{varcovcubpq}}, \code{\link{loglikcubpq}}, \code{\link{inibestgama}}, \code{\link{CUB}}
# @references
# Piccolo D. and D'Elia A. (2008), A new approach for modelling consumers' preferences, \emph{Food Quality and Preference},
# \bold{18}, 247--259 \cr
# Iannario M. and Piccolo D. (2010), A new statistical model for the analysis of customer satisfaction, 
# \emph{Quality Technology and Quantitative Management}, \bold{17}(2)  149--168
#' @keywords internal 

cubpq <-
function(m,ordinal,Y,W,maxiter,toler,makeplot){
  tt0<-proc.time()
  n<-length(ordinal)
  p<-NCOL(Y)
  q<-NCOL(W)
  aver<-mean(ordinal)
  YY<-cbind(1,Y);     WW<-cbind(1,W);          
  #################################################################################
  freq<-tabulate(ordinal,nbins=m)
  inipaicsi<-inibest(m,freq); pai<-inipaicsi[1]; bet0<-log(pai/(1-pai)); betjj<-c(bet0,rep(0.1,p));
  gamajj<-inibestgama(m,ordinal,W)  ### Attention !!!
  #################################################################################
  loglikjj<-loglikcubpq(m,ordinal,Y,W,betjj,gamajj)
  # ********************************************************************
  # ************* E-M algorithm for CUB(p,q) ***************************
  # ********************************************************************
  nniter<-1
  while(nniter<=maxiter){
    loglikold<-loglikjj
    vettn<-bitgama(m,ordinal,W,gamajj)   
    aai<- -1+1/(logis(Y,betjj))   
    ttau<-1/(1+aai/(m*vettn))        
    ####################  maximize w.r.t. bet and gama    ############
    esterno10<-cbind(ttau,YY)
    esterno01<-cbind(ttau,ordinal,WW)
    bet<-betjj;  gama<-gamajj;
    betoptim<-optim(bet,effe10,esterno10=esterno10)
    gamaoptim<-optim(gama,effe01,esterno01=esterno01,m=m)
    ################################################################         
    betjj<-betoptim$par
    gamajj<-gamaoptim$par
    loglikjj<-loglikcubpq(m,ordinal,Y,W,betjj,gamajj)
    # print(c(nniter,betjj,gamajj,loglikjj)); #OPTIONAL PRINTING OF ITERATIONS
    testll<-abs(loglikjj-loglikold)
    if(testll<=toler) break else {loglikold<-loglikjj}
    nniter<-nniter+1
  }
  bet<-betjj;  gama<-gamajj;  loglik<-loglikjj;
  ####################################################################
  AICCUBpq<- -2*loglik+2*(p+q+2)
  BICCUBpq<- -2*loglik+log(n)*(p+q+2)
  ####################################################################
  # Compute asymptotic standard errors of ML estimates
  ####################################################################
  varmat<-varcovcubpq(m,ordinal,Y,W,bet,gama)
  #if(det(varmat)<=0) stop("Variance-covariance matrix NOT positive definite")
  nomi<-c(paste("beta",0:(length(bet)-1),sep="_"),paste("gamma",0:(length(gama)-1),sep="_"))
  stime<-round(c(bet,gama),5)
  nparam<-length(stime)
  if (isTRUE(varmat==matrix(NA,nrow=nparam,ncol=nparam))==TRUE){
    ddd<-cormat<-matrix(NA,nrow=nparam,ncol=nparam)
    ICOMP<-trvarmat<-NA
    errstd<-wald<-pval<-rep(NA,nparam)
  } else {
    ddd<-diag(sqrt(1/diag(varmat)))
    cormat<-(ddd%*%varmat)%*%ddd  
    trvarmat<-sum(diag(varmat))
    ICOMP<- -2*loglik + nparam*log(trvarmat/nparam) - log(det(varmat))
    errstd<-round(sqrt(diag(varmat)),5);  wald<-round(stime/errstd,5);
    pval<-round(2*(1-pnorm(abs(wald))),20)
  }
  ####################################################################
  # Print CUB(p,q) results of ML estimation  
  ####################################################################
  cat("\n")
  cat("==========================================================================","\n")
  cat("============== CUB Program: version 4.0 (September 2015) =================","\n")
  cat("==========================================================================","\n")
  cat("=======>>> C U B (p,q) model <<<=========   ML-estimates via E-M algorithm   ","\n")
  cat("==========================================================================","\n")
  cat(" Covariates for pai ==> p=",p,"     and     Covariates for csi ==> q=", q,"\n")
  cat("==========================================================================","\n")
  cat("*** m=", m,"  *** Sample size: n=", n,"   *** Iterations=",nniter,"Maxiter=",maxiter,"\n")
  cat("==========================================================================","\n")
  cat("parameters  ML-estimates  stand.errors    Wald-test      p-value ","\n")
  cat("==========================================================================","\n")
  for(i in 1:length(nomi)){
    cat(nomi[i],"     ",stime[i],"      ",errstd[i],"       ",wald[i],"      ",pval[i],"\n")
  }
  ####################################################################
  cat("==========================================================================","\n")
  cat("                         Parameters correlation matrix","\n") 
  rownames(cormat)<-nomi;colnames(cormat)<-nomi; 
  print(round(cormat,5))
  ##############################################################################
  cat("==========================================================================","\n")
  cat("Log-lik(beta^,gamma^)=",round(loglik,digits=8),"\n")
  cat("Mean Log-likelihood  =",round(loglik/n,digits=8),"\n")
  cat("--------------------------------------------------------------------------","\n")
  cat("AIC-CUBpq          =",round(AICCUBpq,digits=8),"\n")
  cat("BIC-CUBpq          =",round(BICCUBpq,digits=8),"\n")
  cat("ICOMP-CUBpq        =",round(ICOMP,digits=8),"\n")
  cat("==========================================================================","\n")
  ### Comparing and plotting distributions if the same covariate for pai and csi 
  ### is dichotomus (0,1) #####
  if(p==1 & length(unique(Y))==2 & q==1 & length(unique(W))==2){
    #### code for dicopaicsi  
    vett<-as.matrix(c(0,1))
    pai0<-logis(vett[1],bet)
    pai1<-logis(vett[2],bet)
    csi0<-logis(vett[1],gama)
    csi1<-logis(vett[2],gama)
    prob0<-probcub00(m,pai0,csi0)
    prob1<-probcub00(m,pai1,csi1)
    maxpr<-max(prob0,prob1)
    if(makeplot==TRUE){
      plot(1:m,prob0,ylim=c(0.0,1.1*maxpr),cex.main=0.8,las=1,
           main="CUB distributions, given pai and csi covariates=0, 1",cex=2,
           xlab="",ylab="Prob(R|D=0)  and  Prob(R|D=1)",pch=1,lty=1,type="b");
      lines(1:m,prob1,cex=2,pch=19,lty=2,type="b");
      abline(h=0);
    }
    ### Expected moments given D=0,1 
    exp0<-expcub00(m,pai0,csi0);      exp1<-expcub00(m,pai1,csi1);
    cubmode0<-which.max(prob0);              cubmode1<-which.max(prob1);
    ### Sample averages and modal value, given D=0,1
    ord0<-ordinal[Y==0]; ord1<-ordinal[Y==1];
    n0<-length(ord0);    n1<-length(ord1);
    aver0<-mean(ord0);   aver1<-mean(ord1);
    obsmode0<-which.max(tabulate(ordinal[Y==0]))
    obsmode1<-which.max(tabulate(ordinal[Y==1]))
    cat("Samples and populations measures, given dichotomous covariate (D=0) and (D=1)","\n") 
    cat("----------------------------------------------------------------------------","\n")
    cat("(D = 0)","   n0 = ", n0,
        "        pai_0=",round(pai0,digits=3),"   csi_0=",round(csi0,digits=3),"\n")
    cat("............................","\n")
    cat("Sample average  =",round(aver0,digits=8),"   Sample mode =",round(obsmode0,digits=1),"\n")
    cat("CUB expectation =",round(exp0,digits=8), "    CUB mode   =",round(cubmode0,digits=1),"\n")
    cat("----------------------------------------------------------------------------","\n")
    cat("(D = 1)","   n1 = ", n1,
        "        pai_1=",round(pai1,digits=3),"   csi_1=",round(csi1,digits=3),"\n")
    cat("............................","\n")
    cat("Sample average  =",round(aver1,digits=8)," Sample mode =",round(obsmode1,digits=1),"\n")
    cat("CUB expectation =",round(exp1,digits=8), " CUB mode    =",round(cubmode1,digits=1),"\n")
    cat("----------------------------------------------------------------------------","\n")
  }
  ################################################################
  #        Assignments as global variables
  ################################################################
  #   assign('bet',bet,pos=1)
  #   assign('gama',gama,pos=1)
  #   assign('varmat',varmat,pos=1)
  #   assign('loglik',loglik,pos=1)
  #   assign('niter',nniter,pos=1)
  ################################################################
  durata<-proc.time()-tt0;durata<-durata[1];
  cat("==========================================================================","\n")
  cat("Convergence code for gamma=",gamaoptim$convergence,"\n")
  cat("Convergence code for beta =",betoptim$convergence,"\n")
  cat("==========================================================================","\n")
  cat("Elapsed time              =",durata,"seconds","=====>>>",date(),"\n")
  results<-list('estimates'=stime,'loglik'=loglik,'niter'=nniter,'varmat'=varmat, 'BIC'=round(BICCUBpq,digits=8))
}
