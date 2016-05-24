
pnpmle =function(n,t=15,C=0,b=200,seed=NULL,conf=0.95,dis=1){

  ## ==================================================================================
  ## Purpose: This function calculates the penalized NPMLE by Wang and Lindsay 2005
  ## input:   n--- frequency and frequency data n
  ##          t --- integer cutoff value defining less abundant species, default is 15
  ##          conf--confidence level, a numerical value<1, default .95
  ##          seed--random seed for bootstrap.
  ##          bootrap--number of bootstrap samples.
  ##          dis--0 or 1, 0 for NO display on screen, 1 for yes.
  ## output:  Point estimator and bootstrap confidence interval.
  ## =================================================================================


  if (t!=round(t)||t<0) stop("Error: The cutoff t to define less abundant species must be non-negative integer!")
  if (ncol(n)!=2||is.numeric(n[,1])==FALSE||is.numeric(n[,2])==FALSE) {
    stop("Error: The frequency of frequencies data n must be a matrix.")
  }
  if (C!=0 && C!=1){
    stop("Error: The argument C must be equal to 0 (w/o confidence interval output) or 1 (with confidence inerval)")
  }
  if(is.null(seed)==FALSE) set.seed(seed)
  if(is.numeric(conf)==FALSE||conf>1||conf<0) stop("Error: confidence level must be a numerical value between 0 and 1, e.g. 0.95")

  
  n=as.matrix(n)
  m=max(n[,1])
  ntemp=cbind(c(1:m),rep(0,m))
  ntemp[n[,1],2]=n[,2]
  n=as.matrix(ntemp)
  colnames(n)=NULL; rownames(n)=NULL 
  CI0=numeric(2);lb=0;ub=0;MLE=-1.0
   
  if(t<=nrow(n) && nrow(n)<=50){	
    ntemp=c(n[1:t,2],rep(0,50-t))
    while(MLE<0){
      MLE=-1.0
      p=rep(0.,10)
      pi=p
      noZeroP=0
      theta0=0.0
      results=.Fortran("Nwl",as.double(ntemp),as.integer(t),MLE=as.double(MLE),p=as.double(p),pi=as.double(pi),noZeroP=as.integer(noZeroP),PACKAGE="SPECIES")
      MLE=results$MLE
    }
    MLE0=MLE+sum(n[,2])-sum(n[1:t,2])
    p=results$p
    pi=results$pi
    noZeroP=results$noZeroP
    theta0=results$theta0

    temp=list(MLE0=MLE0,p=p,pi=pi,noZeroP=noZeroP,dis=dis)
    class(temp)="pnpmleClass"
    print(temp)

    ##bootstrap CI
    
    if(C==1){
      rep=as.integer(MLE)
      results=bootPnpmle(p[1:noZeroP],pi[1:noZeroP],b,t,rep,conf)
      lb=results[[1]]+sum(n[,2])-sum(n[1:t,2])
      ub=results[[2]]+sum(n[,2])-sum(n[1:t,2])
    }
  } else if (t>nrow(n)&&nrow(n)<=50){
    t=nrow(n)
    ntemp=c(n[1:t,2],rep(0,50-t))
    while(MLE<0){
      MLE=-1.0
      p=rep(0.,10)
      pi=p
      noZeroP=0
      theta0=0.0
      results=.Fortran("Nwl",as.double(ntemp),as.integer(t),MLE=as.double(MLE),p=as.double(p),pi=as.double(pi),noZeroP=as.integer(noZeroP),PACKAGE="SPECIES")
      MLE=results$MLE
    }
    
    MLE=results$MLE
    MLE0=MLE+sum(n[,2])-sum(n[1:t,2])
    p=results$p
    pi=results$pi
    noZeroP=results$noZeroP
    theta0=results$theta0

    temp=list(MLE0=MLE0,p=p,pi=pi,noZeroP=noZeroP,dis=dis)
    class(temp)="pnpmleClass"
    print(temp)
 
    ##bootstrap CI
    if(C==1){
      rep=as.integer(MLE)
      results=bootPnpmle(p[1:noZeroP],pi[1:noZeroP],b,t,rep,conf)
      lb=results[[1]]
      ub=results[[2]]
    }
  } else if (nrow(n)>50){
    if(t>50){
      ##define less abundant species using t=50
      t=50
    }
    
    ntemp=c(n[1:t,2])
    while(MLE<0){
      MLE=-1.0
      p=rep(0.,10)
      pi=p
      noZeroP=0
      theta0=0.0
      results=.Fortran("Nwl",as.double(ntemp),as.integer(t),MLE=as.double(MLE),p=as.double(p),pi=as.double(pi),noZeroP=as.integer(noZeroP),PACKAGE="SPECIES")
      MLE=results$MLE
    }  
    MLE0=MLE+sum(n[,2])-sum(n[1:t,2])
    p=results$p
    pi=results$pi
    noZeroP=results$noZeroP
    amodel=results$amodel
    amodel0=amodel
    
    temp=list(MLE0=MLE0,p=p,pi=pi,noZeroP=noZeroP,dis=dis)
    class(temp)="pnpmleClass"
    print(temp)

    
    if(C==1){
      rep=as.integer(MLE)
      results=bootPnpmle(p[1:noZeroP],pi[1:noZeroP],b,t,rep,conf)
      lb=results[[1]]+sum(n[,2])-sum(n[1:t,2])
      ub=results[[2]]+sum(n[,2])-sum(n[1:t,2])
    }
  }
  
   if(C==1){
     CI0=matrix(round(c(lb,ub)),1,2)
     colnames(CI0)=c("lb","ub")
     cat("\n")
     return(list(Nhat=round(MLE0),CI=CI0))
   } else if(C==0){
     return(list(Nhat=round(MLE0)))
   }
}


##function to boot pNPMLE
bootPnpmle=function(p,pi,b,t,rep,conf){
  pistat=pi*(1-exp(-p))^(-1)/sum(pi*(1-exp(-p))^(-1))
  SampleProb=untrunPmix(c(0:t),p,pistat)
  SampleProb=SampleProb/(sum(SampleProb))
  Nest=numeric(b)
  cat("\n","Start bootstrap ", b," times:","\n")
  for (i in 1:b){
    MLE=-1.0
    x=sample(c(0:t),rep,replace=TRUE,prob=SampleProb)
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
      results=.Fortran("Nwl",as.double(boot),as.integer(t),MLE=as.double(MLE),p=as.double(p),pi=as.double(pi),noZeroP=as.integer(noZeroP),PACKAGE="SPECIES")
      MLE=results$MLE
    }
    ##cat(boot,"\n")
    Nest[i]=results$MLE
    if(i/20==round(i/20)){
      cat("*")
      if(i/100==round(i/100)){cat("\n")}
    }else{cat(".")}
    
    flush.console()
  }
  cat("\n")
  Nest=sort(Nest)

  lb=round(quantile(Nest,(1-conf)/2))
  ub=round(quantile(Nest,1-(1-conf)/2))
  return(list(lb,ub))
}

#S3 method for on-screen plot
print.pnpmleClass=function(results){ #results is  list of (MLE0,p,pi,noZeroP,dis)
  if(results$dis==1){
    cat("Method: Penalized NPMLE method by  Wang and Lindsay 2005.", "\n\n")
    cat('        MLE=                               ', results$MLE0,"\n")
    cat('        Estimated zero-truncated Poisson mixture components:      ', "\n")
    cat('        p=                                 ', results$p[1:results$noZeroP],"\n")
    cat('        pi=                                ', results$pi[1:results$noZeroP],"\n\n")
  }
}
