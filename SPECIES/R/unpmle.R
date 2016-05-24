
unpmle =function(n,t=15,C=0,method="W-L",b=200,conf=.95,seed=NULL,dis=1){
  
  ## ================================================================================================
  ## Purpose: This function calculates the unconditional NPMLE by Norris and Pollock 1998 using
  ##          algorithm by Bonhing and Schon 2005; or calculates the approximate unconditional NPMLE
  ##          using penalized NPMLE by Wang and Lindsay 2005.
  ## input:   n--- frequency and frequency data n
  ##          t --- integer cutoff value defining less abundant species, default is 15
  ##          method --- string, either method "W-L" or "N-P".
  ##          conf--confidence level, a numerical value<1, default .95
  ##          seed--random seed for bootstrap.
  ##          bootrap--number of bootstrap samples.
  ##          dis--0 or 1, 0 for NO display on screen, 1 for yes.
  ## output:  Point estimator and bootstrap confidence interval.
  ## =================================================================================================
  
  
  if (t!=round(t)||t<0) stop("Error: The cutoff t to define less abundant species must be a non-negative integer!")
  if (ncol(n)!=2||is.numeric(n[,1])==FALSE||is.numeric(n[,2])==FALSE) {
    stop("Error: The frequency of frequencies data n must be a matrix.")
  }
  if (C!=0 && C!=1){
    stop("Error: The argument C must be equal to 0 (w/o confidence interval output) or 1 (with confidence inerval)")
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
  lb=0;ub=0;CI0=numeric(2)
  MLE=-1.0
  
  
  ##cat("N-estimate will be calculated using less abundant species with frequency <= ",t,".","\n\n")  
  
  if(t<=nrow(n) && nrow(n)<=50){	
    ntemp=c(n[1:t,2],rep(0,50-t))*1.0
    if(method=="N-P"){
      while(MLE<0){
        MLE=-1.0
        p=rep(0.,10)
        pi=p
        noZeroP=0
        noZeroP=as.integer(noZeroP)
        results=.Fortran("norrispollock",as.double(ntemp),as.integer(t),MLE=as.double(MLE),p=as.double(p),pi=as.double(pi),noZeroP=as.integer(noZeroP),PACKAGE="SPECIES")
        MLE=results$MLE
      }
      p=results$p
      pi=results$pi
      noZeroP=results$noZeroP
    } else if (method=="W-L"){
      while(MLE<0){
        MLE=-1.0
        p=rep(0.,10)
        pi=p
        noZeroP=0
        noZeroP=as.integer(noZeroP)
        results=.Fortran("WLunpmle",as.double(ntemp),as.integer(t),MLE=as.double(MLE),p=as.double(p),pi=as.double(pi),noZeroP=as.integer(noZeroP),PACKAGE="SPECIES")
        MLE=results$MLE
      }
      p=results$p
      pi=results$pi
      noZeroP=results$noZeroP
      pistat=pi[1:noZeroP]*(1-exp(-p[1:noZeroP]))^(-1)/sum(pi[1:noZeroP]*((1-exp(-p[1:noZeroP]))^(-1)))
      pi[1:noZeroP]=pistat  
    }
    
    MLE0=MLE+sum(n[,2])-sum(n[1:t,2])

    temp=list(MLE0=MLE0,p=p,pi=pi,noZeroP=noZeroP,dis=dis,method=method)
    class(temp)="unpmleClass"
    print(temp)
         
    if(C==1){
      rep=sum(ntemp)
      results=bootUnpmle(p[1:noZeroP],pi[1:noZeroP],b,t,method,rep,conf)
      lb=results[[1]]+sum(n[,2])-sum(n[1:t,2])
      ub=results[[2]]+sum(n[,2])-sum(n[1:t,2])
    }
  } else if (t>nrow(n)&&nrow(n)<=50){
    t=nrow(n)
    ntemp=c(n[1:t,2],rep(0,50-t))*1.0
    if(method=="N-P"){
      while(MLE<0){
        MLE=-1.0
        p=rep(0.,10)
        pi=p
        noZeroP=0
        noZeroP=as.integer(noZeroP)
        results=.Fortran("norrispollock",as.double(ntemp),as.integer(t),MLE=as.double(MLE),p=as.double(p),pi=as.double(pi),noZeroP=as.integer(noZeroP),PACKAGE="SPECIES")
        MLE=results$MLE
      }
      p=results$p
      pi=results$pi
      noZeroP=results$noZeroP
    } else if (method=="W-L"){
      while(MLE<0){
        MLE=-1.0
        p=rep(0.,10)
        pi=p
        noZeroP=0
        results=.Fortran("WLunpmle",as.double(ntemp),as.integer(t),MLE=as.double(MLE),p=as.double(p),pi=as.double(pi),noZeroP=as.integer(noZeroP),PACKAGE="SPECIES")
        
        MLE=results$MLE
      }
      p=results$p
      pi=results$pi
      noZeroP=results$noZeroP
      pistat=pi[1:noZeroP]*(1-exp(-p[1:noZeroP]))^(-1)/sum(pi[1:noZeroP]*((1-exp(-p[1:noZeroP]))^(-1)))
      pi[1:noZeroP]=pistat  
    }
    
    MLE0=MLE+sum(n[,2])-sum(n[1:t,2])
    temp=list(MLE0=MLE0,p=p,pi=pi,noZeroP=noZeroP,dis=dis,method=method)
    class(temp)="unpmleClass"
    print(temp)
    
    if(C==1){
      rep=sum(ntemp)
      results=bootUnpmle(p[1:noZeroP],pi[1:noZeroP],b,t,method,rep,conf)
      lb=results[[1]]+sum(n[,2])-sum(n[1:t,2])
      ub=results[[2]]+sum(n[,2])-sum(n[1:t,2])
    }
  } else if (nrow(n)>50){
    if(t>50){
      ##define less abundant species using t=50
      t=50
    }
    ntemp=c(n[1:t,2])*1.0
    if(method=="N-P"){
      while(MLE<0){
        MLE=-1.0
        p=rep(0.,10)
        pi=p
        noZeroP=0
        noZeroP=as.integer(noZeroP)
        results=.Fortran("norrispollock",as.double(ntemp),as.integer(t),MLE=as.double(MLE),p=as.double(p),pi=as.double(pi),noZeroP=as.integer(noZeroP),PACKAGE="SPECIES")
        MLE=results$MLE
      }
      p=results$p
      pi=results$pi
      noZeroP=results$noZeroP
    } else if (method=="W-L"){
      while(MLE<0){
        MLE=-1.0
        p=rep(0.0,10)
        pi=p
        noZeroP=0
        noZeroP=as.integer(noZeroP)
        results=.Fortran("WLunpmle",as.double(ntemp),as.integer(t),MLE=as.double(MLE),p=as.double(p),pi=as.double(pi),noZeroP=as.integer(noZeroP),PACKAGE="SPECIES")
        
        MLE=results$MLE
      }
      p=results$p
      pi=results$pi
      noZeroP=results$noZeroP
      pistat=pi[1:noZeroP]*(1-exp(-p[1:noZeroP]))^(-1)/sum(pi[1:noZeroP]*((1-exp(-p[1:noZeroP]))^(-1)))
      pi[1:noZeroP]=pistat  
    }
    
    MLE0=MLE+sum(n[,2])-sum(n[1:t,2])
    temp=list(MLE0=MLE0,p=p,pi=pi,noZeroP=noZeroP,dis=dis,method=method)
    class(temp)="unpmleClass"
    print(temp)
    
    if(C==1){
      rep=sum(ntemp)
      results=bootUnpmle(p[1:noZeroP],pi[1:noZeroP],b,t,method,rep,conf)
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
  else{}
}


##bootstrap Unpmle
bootUnpmle=function(p,pi,b,t,method,rep,conf){
  SampleProb=untrunPmix(c(1:t),p,pi)
  SampleProb=SampleProb/(sum(SampleProb))
  Nest=numeric(b)
  cat("\n","Start bootstrap", b,  " times:","\n")
  i=1
  while (i <= b){
    MLE=-1.0
    x=sample(c(1:t),rep,replace=TRUE,prob=SampleProb)

    if(method=="N-P"){
      while(MLE<0){
        boot=numeric(50)
        for(j in 1:t){
          boot[j]=sum(x==j)
        }
        boot=as.double(boot)

        MLE=-1.0
        MLE=as.double(MLE)
        p=rep(0.,10)
        p=as.double(p)
        pi=p
        noZeroP=0
        noZeroP=as.integer(noZeroP)
        results=.Fortran("norrispollock",as.double(boot),as.integer(t),MLE=as.double(MLE),p=as.double(p),pi=as.double(pi),noZeroP=as.integer(noZeroP),PACKAGE="SPECIES")
        MLE=results$MLE
        #cat("i=",i," ",MLE,"\n")
      }
    } else if (method=="W-L"){
      while(MLE<0){
    	boot=numeric(50)
        for(j in 1:t){
          boot[j]=sum(x==j)
        }
        boot=as.double(boot)
        MLE=-1.0
        MLE=as.double(MLE)
        p=rep(0.,10)
        p=as.double(p)
        pi=p
        noZeroP=0
        noZeroP=as.integer(noZeroP)
        results=.Fortran("WLunpmle",as.double(boot),as.integer(t),MLE=as.double(MLE),p=as.double(p),pi=as.double(pi),noZeroP=as.integer(noZeroP),PACKAGE="SPECIES")
        MLE=results$MLE
        #cat("i=",i," ",MLE,"\n")
      }
    }
    Nest[i]=results$MLE
    if(i/20==round(i/20)){
      cat("*")
      if(i/100==round(i/100)){cat("\n")}
    }else{cat(".")}
    
    i=i+1
    ##cat (boot,"\n")
    flush.console()
  }
  cat("\n")
  lb=round(quantile(Nest,(1-conf)/2))
  ub=round(quantile(Nest,1-(1-conf)/2))
  return(list(lb,ub))
}


#S3 method for on-screen plot
print.unpmleClass=function(results){ #results is  list of (MLE0,p,pi,noZeroP,dis,method)
  if(results$dis==1){
    cat("Method: Unconditional NPMLE method by Norris and Pollock 1996, 1998,","\n")
    if(results$method=="N-P"){
      cat("        using algorithm by Bonhing and Schon 2005:", "\n\n")
    }else if(results$method=="W-L"){
      cat("        using algorithm by Wang and Lindsay 2005:", "\n\n")
    }else {
      stop ("Error: the method must be either N-P or W-L!")
    }
 
    cat('        MLE=                               ', results$MLE0,"\n")
    cat('        Estimated Poisson mixture components:      ', "\n")
    cat('        p=                                 ', results$p[1:results$noZeroP],"\n")
    cat('        pi=                                ', results$pi[1:results$noZeroP],"\n\n")
  }
}
