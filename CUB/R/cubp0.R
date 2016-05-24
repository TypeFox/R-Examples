# @title Main function for CUB models with covariates for the uncertainty component
# @description Estimate and validate a CUB model for given ordinal responses, with covariates for explaining 
# the feeling component via a logistic transform.
# @aliases cubp0
# @usage cubp0(m, ordinal, Y, maxiter, toler, makeplot)
# @param m Number of ordinal categories
# @param ordinal Vector of ordinal responses
# @param Y Matrix of selected covariates for explaining the uncertainty component
# @param maxiter Maximum number of iterations allowed for running the optimization algorithm 
# @param toler Fixed error tolerance for final estimates 
# @param makeplot Logical: if TRUE and if only one dichotomous covariate is included in the model, with levels (0,1),  
# the function returns a graphical plot comparing the distributions of the responses conditioned to the value 
# of the covariate
# @return An object of the class "CUB"
# @import stats graphics
# @references
# Iannario M. and Piccolo D. (2010), A new statistical model for the analysis of customer satisfaction, 
# \emph{Quality Technology and Quantity management}, \bold{7}(2) 149--168 \cr
# Iannario M. and Piccolo D. (2012). CUB models: Statistical methods and empirical evidence, in:
#  Kenett R. S. and Salini S. (eds.), \emph{Modern Analysis of Customer Surveys: with applications using R}, 
#  J. Wiley and Sons, Chichester, 231--258
#' @keywords internal 

cubp0 <-
function(m,ordinal,Y,maxiter,toler,makeplot){
  tt0<-proc.time()
  n<-length(ordinal)
  p<-NCOL(Y)
  aver<-mean(ordinal); varcamp<-mean(ordinal^2)-aver^2;
  YY<-cbind(1,Y)
  ##################################################################
  serie<-1:m; freq<-tabulate(ordinal,nbins=m);
  inipaicsi<-inibest(m,freq)
  pai<-inipaicsi[1]; bet0<-log(pai/(1-pai));  
  betjj<- c(bet0,rep(0.1,p))                #betjj<-rep(0.1,p+1);
  csijj<-inipaicsi[2]
  ##############################################################
  loglikjj<-loglikcubp0(m,ordinal,Y,betjj,csijj)
  # ********************************************************************
  # ************* E-M algorithm for CUB(p,0) ***************************
  # ********************************************************************
  nniter<-1
  while(nniter<=maxiter){
    loglikold<-loglikjj
    bb<-probbit(m,csijj)
    vettn<-bb[ordinal]      # probbit for all ordinal (r_i,i=1,2,...,n)
    aai<- -1+ 1/(logis(Y,betjj)) #exp(-(YY%*%betjj));
    ttau<-1/(1+aai/(m*vettn))       # tau is a reserved word in R
    averpo<-sum(ordinal*ttau)/sum(ttau)
    ################################## maximize w.r.t. bet  ########
    bet<-betjj
    covar<-YY
    tauno<-ttau
    #nlmaxbet<-nlm(effe10,betjj,esterno10);   
    opmaxbet<-optim(bet,effe10,esterno10=cbind(tauno,covar))
    ################################################################         
    betjj<-opmaxbet$par
    # betjj<-nlmaxbet$estimate;        #updated bet estimates
    csijj<-(m-averpo)/(m-1)       #updated csi estimate
    #loglikjj<- -opmaxbet$value
    loglikjj<-loglikcubp0(m,ordinal,Y,betjj,csijj)
    
    #print(c(nniter,betjj,csijj,loglikjj)); #OPTIONAL PRINTING OF ITERATIONS
    testll<-abs(loglikjj-loglikold)
    if(testll<=toler) break else {loglikold<-loglikjj}
    nniter<-nniter+1
  }
  bet<-betjj;  csi<-csijj;  loglik<-loglikjj;
  ####################################################################
  AICCUBp0<- -2*loglik+2*(p+2)
  BICCUBp0<- -2*loglik+log(n)*(p+2)
  ####################################################################
  # Compute asymptotic standard errors of ML estimates
  ####################################################################
  varmat<-varcovcubp0(m,ordinal,Y,bet,csi) 
  nomi<-c(paste("beta",0:(length(bet)-1),sep="_"),"csi   ")
  stime<-round(c(bet,csi),5)
  nparam<-length(stime)
  #if(det(varmat)<=0) stop("Variance-covariance matrix NOT positive definite")
  if (isTRUE(varmat==matrix(NA,nrow=nparam,ncol=nparam))==TRUE){
    ddd<-cormat<-matrix(NA,nrow=nparam,ncol=nparam)
    ICOMP<-trvarmat<-NA
    errstd<-wald<-pval<-rep(NA,nparam) 
  } else {
    ddd<-diag(sqrt(1/diag(varmat)))
    cormat<-(ddd%*%varmat)%*%ddd
    trvarmat<-sum(diag(varmat))
    ICOMP<- -2*loglik + nparam*log(trvarmat/nparam) - log(det(varmat)) ## added
    errstd<-round(sqrt(diag(varmat)),5);  wald<-round(stime/errstd,5);
    pval<-round(2*(1-pnorm(abs(wald))),20)
  }
  
  ####################################################################
  # Print CUB(p,0) results of ML estimation  
  ####################################################################
  cat("\n")
  cat("=======================================================================","\n")
  cat("============== CUB Program: version 4.0 (September 2015) ==============","\n")
  cat("=======================================================================","\n")
  cat("=====>>> C U B (p,0) model <<<=====   ML-estimates via E-M algorithm   ","\n")
  cat("=======================================================================","\n")
  cat("                    Covariates for pai ==> p=", p,"\n")
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
  cat("Log-lik(beta^,csi^) =",round(loglik,digits=8),"\n")
  cat("Mean Log-likelihood =",round(loglik/n,digits=8),"\n")
  cat("-----------------------------------------------------------------------","\n")
  cat("AIC-CUBp0          =",round(AICCUBp0,digits=8),"\n")
  cat("BIC-CUBp0          =",round(BICCUBp0,digits=8),"\n")
  cat("ICOMP-CUBp0        =",round(ICOMP,digits=8),"\n")
  cat("=======================================================================","\n")
  ################################################################
  #        Assignments as global variables
  ################################################################
  #   assign('bet',bet,pos=1)
  #   assign('csi',csi,pos=1)
  #   assign('varmat',varmat,pos=1)
  #   assign('loglik',loglik,pos=1)
  #   assign('nniter',nniter,pos=1)
  ### Comparing and plotting distributions if covariate for pai is dichotomus (0,1) #####
  if(p==1 & length(unique(Y))==2) {
    #code dicopai
    vett<-as.matrix(c(0,1))
    pai0<-logis(vett[1],bet)
    pai1<-logis(vett[2],bet)
    prob0<-probcub00(m,pai0,csi)
    prob1<-probcub00(m,pai1,csi)
    maxpr<-max(prob0,prob1)
    if(makeplot==TRUE){
      plot(1:m,prob0,ylim=c(0.0,1.1*maxpr),cex.main=0.8,las=1,
           main="CUB distributions, given pai-covariate=0, 1",cex=2,
           xlab="",ylab="Prob(R|D=0)  and  Prob(R|D=1)",pch=1,lty=1,type="b");
      lines(1:m,prob1,cex=2,pch=19,lty=2,type="b");
      abline(h=0);
    }
    ### Expected moments given D=0,1 
    exp0<-expcub00(m,pai0,csi);      exp1<-expcub00(m,pai1,csi);
    cubmode0<-which.max(prob0);              cubmode1<-which.max(prob1);
    ### Sample averages and modal value, given D=0,1
    ord0<-ordinal[Y==0]; ord1<-ordinal[Y==1];
    n0<-length(ord0);    n1<-length(ord1);
    aver0<-mean(ord0);   aver1<-mean(ord1);
    obsmode0<-which.max(tabulate(ordinal[Y==0]))    
    obsmode1<-which.max(tabulate(ordinal[Y==1]))
    cat("Samples and populations measures, given dichotomous covariate (D=0) and (D=1)","\n")
    cat("-----------------------------------------------------------------------","\n")
    cat("(D = 0)","   n0 = ", n0,
        "        pai_0=",round(pai0,digits=3),"   csi=",round(csi,digits=3),"\n")
    cat("............................","\n")
    cat("Sample average  =",round(aver0,digits=8),"   Sample mode =",round(obsmode0,digits=1),"\n")
    cat("CUB expectation =",round(exp0,digits=8), "    CUB mode   =",round(cubmode0,digits=1),"\n")
    cat("-----------------------------------------------------------------------","\n")
    cat("(D = 1)","   n1 = ", n1,
        "        pai_1=",round(pai1,digits=3),"   csi=",round(csi,digits=3),"\n")
    cat("............................","\n")
    cat("Sample average  =",round(aver1,digits=8)," Sample mode =",round(obsmode1,digits=1),"\n")
    cat("CUB expectation =",round(exp1,digits=8), " CUB mode    =",round(cubmode1,digits=1),"\n")
    cat("-----------------------------------------------------------------------","\n")
  }
  ################################################################
  durata<-proc.time()-tt0;durata<-durata[1];
  cat("=======================================================================","\n")
  cat("Convergence code  =",opmaxbet$convergence,"\n")
  cat("=======================================================================","\n")
  cat("Elapsed time      =",durata,"seconds","=====>>>",date(),"\n")
  results<-list('estimates'=stime, 'loglik'=loglik,'niter'=nniter,'varmat'=varmat,'BIC'=round(BICCUBp0,digits=8))
}
