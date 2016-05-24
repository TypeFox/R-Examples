
pcg =function(n,t=35,C=0,alpha=c(1:10),b=200,seed=NULL,conf=0.95,dis=1){

  ## ===================================================================================
  ## Purpose: This function calculates the Poisson-compound Gamma estimator by Wang 2010
  ## input:   n--- frequency and frequency data n
  ##          t --- integer cutoff value defining less abundant species, default is 35
  ##          conf--confidence level, a numerical value<1, default .95
  ##          seed--random seed for bootstrap.
  ##          bootrap--number of bootstrap samples.
  ##          dis--0 or 1, 0 for NO display on screen, 1 for yes.
  ## output:  Point estimator and bootstrap confidence interval.
  ## ===================================================================================
  
  if (t!=round(t)||t<0) stop("Error: The cutoff t to define less abundant species must be non-negative integer!")
  if (ncol(n)!=2||is.numeric(n[,1])==FALSE||is.numeric(n[,2])==FALSE) {
    stop("Error: The frequency of frequencies data n must be a matrix.")
  }
  if (C!=0 && C!=1){
    stop("Error: The argument C must be equal to 0 (w/o confidence interval output) or 1 (with confidence inerval)")
  }
  if(is.numeric(alpha)==FALSE||sum(alpha<0)>0){
    stop("Error: The argument alpha must be a numeric vector specifying grid points of alpha. The default in 1:10.")
  }
  if(max(alpha)>=500){
    stop("Error: The largest alpha value should be <500.")
  }
  if(is.numeric(conf)==FALSE||conf>1||conf<0) stop("Error: confidence level must be a numerical value between 0 and 1, e.g. 0.95")
  if(is.null(seed)==FALSE) set.seed(seed)
  
  flush.console()
  n=as.matrix(n)
  m=max(n[,1])
  ntemp=cbind(c(1:m),rep(0,m))
  ntemp[n[,1],2]=n[,2]
  n=as.matrix(ntemp)
  colnames(n)=NULL;rownames(n)=NULL
  alphaK=length(alpha)
  
  CI0=numeric(2);lb=0;ub=0; amodel=0.;
  MLE=-1.0
  
  if(t<=nrow(n) && nrow(n)<=50){	
    ntemp=c(n[1:t,2],rep(0,50-t))
    while(MLE<0){
      MLE=-1.0
      p=rep(0.,10)
      pi=p
      noZeroP=0
      noZeroP=as.integer(noZeroP)
      results=.Fortran("PPCG",as.double(ntemp),as.integer(t),MLE=as.double(MLE),as.double(alpha),as.integer(alphaK),amodel=as.double(amodel),p=as.double(p),pi=as.double(pi),noZeroP=as.integer(noZeroP),PACKAGE="SPECIES")
      MLE=results$MLE
    }
    MLE0=MLE+sum(n[,2])-sum(n[1:t,2])
    p=results$p
    pi=results$pi
    noZeroP=results$noZeroP
    amodel=results$amodel
    amodel0=amodel
    if (amodel0==500) amodel0=Inf
    temp=list(MLE0=MLE0,alpha=alpha,amodel0=amodel0,p=p,pi=pi,noZeroP=noZeroP,dis=dis)
    class(temp)="pcgClass"
    print(temp)
   
    if(C==1){
      rep=as.integer(MLE)
      results=bootPCG(p[1:noZeroP],pi[1:noZeroP],b,t,rep,amodel0,alpha,conf)
      lb=results[[1]]+sum(n[,2])-sum(n[1:t,2])
      ub=results[[2]]+sum(n[,2])-sum(n[1:t,2])
    }
  }else if (t>nrow(n)&&nrow(n)<=50){
    t=nrow(n)
    
    ntemp=c(n[1:t,2],rep(0,50-t))
    while(MLE<0){
      MLE=-1.0
      p=rep(0.,10)
      pi=p
      noZeroP=0
      noZeroP=as.integer(noZeroP)
      results=.Fortran("PPCG",as.double(ntemp),as.integer(t),MLE=as.double(MLE),as.double(alpha),as.integer(alphaK),amodel=as.double(amodel),p=as.double(p),pi=as.double(pi),noZeroP=as.integer(noZeroP),PACKAGE="SPECIES")
      MLE=results$MLE
    }
    
    MLE0=MLE+sum(n[,2])-sum(n[1:t,2])
    p=results$p
    pi=results$pi
    noZeroP=results$noZeroP
    amodel=results$amodel
    amodel0=amodel
    if (amodel0==500) amodel0=Inf
    temp=list(MLE0=MLE0,alpha=alpha,amodel0=amodel0,p=p,pi=pi,noZeroP=noZeroP,dis=dis)
    class(temp)="pcgClass"
    print(temp)
    
    if(C==1){
      rep=sum(n[1:t,2])
      results=bootPCG(p[1:noZeroP],pi[1:noZeroP],b,t,rep,amodel0,alpha,conf)
      lb=results[[1]]+sum(n[,2])-sum(n[1:t,2])
      ub=results[[2]]+sum(n[,2])-sum(n[1:t,2])
    }
  } else if (nrow(n)>50){
    if(t>50){
      ##define less abundant species using t=50
      t=50
    }
    ntemp=c(n[1:t,2])
    MLE=-1.0
    while(MLE<0){
      MLE=-1.0
      p=rep(0.,10)
      pi=p
      noZeroP=0
      noZeroP=as.integer(noZeroP)
      results=.Fortran("PPCG",as.double(ntemp),as.integer(t),MLE=as.double(MLE),as.double(alpha),as.integer(alphaK),amodel=as.double(amodel),p=as.double(p),pi=as.double(pi),noZeroP=as.integer(noZeroP),PACKAGE="SPECIES")
      MLE=results$MLE
    }
    MLE0=MLE+sum(n[,2])-sum(n[1:t,2])
    p=results$p
    pi=results$pi
    noZeroP=results$noZeroP
    amodel=results$amodel
    amodel0=amodel
    if (amodel0==500) amodel0=Inf
    temp=list(MLE0=MLE0,alpha=alpha,amodel0=amodel0,p=p,pi=pi,noZeroP=noZeroP,dis=dis)
    class(temp)="pcgClass"
    print(temp)
    
    if(C==1){
      rep=sum(n[1:t,2])
      results=bootPCG(p[1:noZeroP],pi[1:noZeroP],b,t,rep,amodel0,alpha,conf)
      lb=results[[1]]+sum(n[,2])-sum(n[1:t,2])
      ub=results[[2]]+sum(n[,2])-sum(n[1:t,2])
    }
  }
  
  if(C==1){
    cat("\n")
    CI0=matrix(round(c(lb,ub)),1,2)
    colnames(CI0)=c("lb","ub")
    cat("\n")
    return(list(Nhat=round(MLE0),AlphaModel=amodel0,CI=CI0))
  }
  else if(C==0){
    return(list(Nhat=round(MLE0),AlphaModel=amodel0))
  }
  else{}
}

  

###
bootPCG= function(p,pi,b,t,rep,amodel0,alpha,conf){
  
  if(amodel0<Inf){
    SampleProb=PmixSCon(c(1:t),p,pi,amodel0)
  } else{
    SampleProb=Pmix(c(1:t),p,pi)
  }
  Nest=numeric(b)
  alphaK=length(alpha)
  amodel=0
  
  cat("\n","Start bootstrap", b,  " times:","\n")
  
  for (i in 1:b){
    ##cat(boot)
    MLE=-1.0
    x=sample(c(1:t),rep,replace=TRUE,prob=SampleProb)

    while(MLE<0){
      boot=numeric(50)
      for(j in 1:t){
        boot[j]=sum(x==j)
      }
      
      MLE=-1.0
      p=rep(0.,10)
      pi=p
      noZeroP=0
      theta0=0.0
      results=.Fortran("PPCG",as.double(boot),as.integer(t),MLE=as.double(MLE),as.double(alpha),as.integer(alphaK),amodel=as.double(amodel),p=as.double(p),pi=as.double(pi),noZeroP=as.integer(noZeroP),PACKAGE="SPECIES")
      MLE=results$MLE
    }
    Nest[i]=results$MLE
    if(i/20==round(i/20)){
      cat("*")
      if(i/100==round(i/100)){cat("\n")}
    }else{cat(".")}
    flush.console()
  }
  cat("\n")

  lb=round(quantile(Nest,(1-conf)/2))
  ub=round(quantile(Nest,1-(1-conf)/2))
  return(list(lb,ub))
}

###s3 print method
print.pcgClass=function(results){#results is  list of (MLE0,amodel0,alpha,p,pi,noZeroP,dis)
  if(results$dis==1){
    cat("Method: Poisson-Compound Gamma method by Wang 2010.", "\n")
    cat("Alpha grid used:", results$alpha,".","\n\n")
    cat('        MLE=                               ', results$MLE0,"\n")
    cat('        Selected alpha model:              ', results$amodel0,"\n")
    cat('        Estimated Gamma components:        ', "\n")
    cat('        p=                                 ', results$p[1:results$noZeroP],"\n")
    cat('        pi=                                ', results$pi[1:results$noZeroP],"\n\n")
  }
}
