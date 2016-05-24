
#' @keywords internal 


cub0q <-function(m,ordinal,W,maxiter,toler,makeplot){
  tt0<-proc.time()
  n<-length(ordinal)
  q<-NCOL(W)
  aver<-mean(ordinal)
  WW<-cbind(1,W)                   
  ##############################################################
  freq<-tabulate(ordinal,nbins=m); inipaicsi<-inibest(m,freq);  paijj<-inipaicsi[1]; 
  gamajj<-inibestgama(m,ordinal,W)
  ##############################################################
  loglikjj<-loglikcub0q(m,ordinal,W,paijj,gamajj)
  # ********************************************************************
  # ************* E-M algorithm for CUB(0,q) ***************************
  # ********************************************************************
  nniter<-1
  while(nniter<=maxiter){
    loglikold<-loglikjj
    vettn<-bitgama(m,ordinal,W,gamajj)
    ttau<-1/(1+(1-paijj)/(m*paijj*vettn)) 
    ################################# maximize w.r.t. gama  ########
    ordd<-ordinal;covar<-WW;
    gama<-gamajj
    optimgama<-optim(gama,effe01,esterno01=cbind(ttau,ordinal,WW),m=m) 
    ################################################################
    gamajj<-optimgama$par
    paijj<-sum(ttau)/n                    #updated pai estimate
    loglikjj<-loglikcub0q(m,ordinal,W,paijj,gamajj)## needed for nlm version
    # print(c(nniter,paijj,gamajj,loglikjj)); #OPTIONAL PRINTING OF ITERATIONS
    testll<-abs(loglikjj-loglikold)
    if(testll<=toler) break else {loglikold<-loglikjj}
    nniter<-nniter+1
  }
  pai<-paijj;  gama<-gamajj;  loglik<-loglikjj;
  ####################################################################
  AICCUB0q<- -2*loglik+2*(q+2)
  BICCUB0q<- -2*loglik+log(n)*(q+2)
  ####################################################################
  # Compute asymptotic standard errors of ML estimates
  ####################################################################
  varmat<-varcovcub0q(m,ordinal,W,pai,gama)
  nomi<-c("pai    ",paste("gamma",0:(length(gama)-1),sep="_"))
  stime<-round(c(pai,gama),5)
  nparam<-length(stime)
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
  ### Print CUB(0,q) results of ML estimation  
  ####################################################################
  cat("\n")
  cat("=======================================================================","\n")
  cat("============== CUB Program: version 4.0 (September 2015) ==============","\n")
  cat("=======================================================================","\n")
  cat("=====>>> C U B (0,q) model <<<=====   ML-estimates via E-M algorithm   ","\n")
  cat("=======================================================================","\n")
  cat("                    Covariates for csi ==> q=", q,"\n")
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
  cat("Log-lik(pai^,gamma^) =",round(loglik,digits=8),"\n")
  cat("Mean Log-likelihood  =",round(loglik/n,digits=8),"\n")
  cat("-----------------------------------------------------------------------","\n")
  cat("AIC-CUB0q          =",round(AICCUB0q,digits=8),"\n")
  cat("BIC-CUB0q          =",round(BICCUB0q,digits=8),"\n")
  cat("ICOMP-CUB0q        =",round(ICOMP,digits=8),"\n")
  cat("=======================================================================","\n")
  ################################################################
  #        Assignments as global variables
  ################################################################
  #   assign('pai',pai,pos=1)
  #   assign('gama',gama,pos=1)
  #   assign('varmat',varmat,pos=1)
  #   assign('loglik',loglik,pos=1)
  #   assign('nniter',nniter,pos=1)
  ### Comparing and plotting distributions if covariate for csi is dichotomus (0,1) #####
  if(q==1 & length(unique(W))==2) {
    #code for dicocsi
    vett<-as.matrix(c(0,1))
    csi0<-logis(vett[1],gama); prob0<-probcub00(m,pai,csi0);
    csi1<-logis(vett[2],gama); prob1<-probcub00(m,pai,csi1);
    maxpr<-max(prob0,prob1)
    ## makeplot=TRUE ### Alternatively, makeplot=FALSE
    if(makeplot==TRUE){
      plot(1:m,prob0,ylim=c(0.0,1.1*maxpr),cex.main=0.8,las=1,
           main="CUB distributions, given csi-covariate=0, 1",
           cex=2,xlab="",ylab="Prob(R|D=0)  and  Prob(R|D=1)",pch=1,lty=1,type="b");
      lines(1:m,prob1,cex=2,pch=19,lty=2,type="b");
      abline(h=0);
    }
    ### Expected moments given D=0,1 
    exp0<-expcub00(m,pai,csi0);      exp1<-expcub00(m,pai,csi1);
    cubmode0<-which.max(prob0);      cubmode1<-which.max(prob1);
    ### Sample averages and modal value, given D=0,1
    ord0<-ordinal[W==0]; ord1<-ordinal[W==1];
    n0<-length(ord0);    n1<-length(ord1);
    aver0<-mean(ord0);   aver1<-mean(ord1);
    obsmode0<-which.max(tabulate(ord0))
    obsmode1<-which.max(tabulate(ord1))
    cat("Samples and populations measures, given dichotomous covariate (D=0) and (D=1)","\n")
    cat("-----------------------------------------------------------------------","\n")
    cat("(D = 0)","   n0 = ", n0,
        "        pai=",round(pai,digits=3),"   csi_0=",round(csi0,digits=3),"\n")
    cat("............................","\n")
    cat("Sample average  =",round(aver0,digits=8),"   Sample mode =",round(obsmode0,digits=1),"\n")
    cat("CUB expectation =",round(exp0,digits=8), "   CUB mode    =",round(cubmode0,digits=1),"\n")
    cat("-----------------------------------------------------------------------","\n")
    cat("(D = 1)","   n1 = ", n1,
        "        pai=",round(pai,digits=3),"   csi_1=",round(csi1,digits=3),"\n")
    cat("............................","\n")
    cat("Sample average  =",round(aver1,digits=8)," Sample mode =",round(obsmode1,digits=1),"\n")
    cat("CUB expectation =",round(exp1,digits=8), " CUB mode    =",round(cubmode1,digits=1),"\n")
    cat("-----------------------------------------------------------------------","\n")
  }
  ################################################################
  durata<-proc.time()-tt0;durata<-durata[1];
  cat("=======================================================================","\n")
  cat("Convergence code =",optimgama$convergence,"\n")
  cat("=======================================================================","\n")
  cat("Elapsed time     =",durata,"seconds","=====>>>",date(),"\n")
  results<-list('estimates'=stime,'loglik'=loglik,'niter'=nniter,'varmat'=varmat,'BIC'=round(BICCUB0q,digits=8))
}
