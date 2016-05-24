simSeerSet<-function(N=2e9,yearEnd=2012,ka=1e-5,kb=0.04,Ab=1e-5,tauA=10,tauB=1,delay=1,period=4) {
  #   agedx=age=age86=canc=yrdx=sex=race=surv=modx=yrbrth=NULL 
  trt=cancers=NULL 
  data(stdUS, envir = environment())
  #   if (shape<1) {
  #     shape=1
  #     print("Warning: shape must be 1 or higher. It was less and is being set to 1. ")
  #   }
  #   N=2e9;yearEnd=2012;ka=1e-5;kb=0.04;Ab=1e-5;tauA=10;tauB=5;delay=1;period=4;shape=1;library(dplyr)
  popsa=merge(data.frame(age=0.5:99.5),data.frame(year=1973:yearEnd))
  popsa$py=N/40*SEERaBomb::stdUS$prop[round(popsa$age+0.5)]
  #   sum(popsa$py)
  head(popsa)
  #   tail(popsa)
  A=cbind(popsa[,1:2],cancers=rpois(dim(popsa)[1],ka*popsa$age*popsa$py))
  head(A)
  sum(A$cancers)
  cancA=A[rep(seq_len(nrow(A)), times=A$cancers),]%>%select(-cancers)
  cancA$surv=rexp(dim(cancA)[1],rate=1/tauA)
  cancA$yrdx=cancA$year+runif(dim(cancA)[1],max=0.9999)
  head(cancA,10)
  #   cancA$surv=runif(dim(cancA)[1],0,20)
  #   mean(cancA$surv)
  trimSurv=function(x) {
    x$surv=ifelse(x$surv+x$age>123,123-x$age,x$surv)
    x$surv=ifelse(x$surv+x$yrdx>yearEnd+0.9999,yearEnd+0.9999-x$yrdx,x$surv)
    x
  }
  cancA=trimSurv(cancA)
  cancA$cancer="A"
  trts=c("noRad","rad")
  cancA$trt=sample(trts,dim(cancA)[1],replace=T)
  cancA$casenum=1:dim(cancA)[1]
  cancA$seqnum=0
  table(cancA$trt)
  #   rownames(cancA)=cancA$casenum  # can't do as rownames since IDs are not unique
  head(cancA)
  
  ###### treatment of A with rad induces B
  cancAr=cancA%>%filter(trt=="rad")
  py=cancAr$surv
  py=ifelse(py<=delay,0,py-delay)
  py=ifelse(py>=period,period,py)
  cancAr$seconds=rpois(dim(cancAr)[1],4*(Ab/kb)*(exp(kb*(cancAr$age+py+delay))-exp(kb*(cancAr$age+delay)) ) )
  sum(cancAr$seconds)
  cancAr$seqnum[cancAr$seconds>0]=1
  cancBAr=cancAr[cancAr$seconds>0,]
  cancBAr$seqnum=2
  cancBAr$cancer="B"
  cancBAr$trt="noRad"
  
  py=py[cancAr$seconds>0]
  Y=runif(dim(cancBAr)[1],Ab*exp(kb*(cancBAr$age+delay)),Ab*exp(kb*(cancBAr$age+py+delay)))
  ages=(log(Y)-log(Ab))/kb
  cancBAr$yrdx=cancBAr$yrdx+ages-cancBAr$age
  
  
  #   cancBA$yrdx=cancBA$yrdx+delay+period*rbeta(dim(cancBA)[1],shape,shape)
  #   cancBA$yrdx=cancBA$yrdx+runif(dim(cancBA)[1],delay+1e-4,delay+period-1e-4)
  #   cancBA$yrdx=cancBA$yrdx+     runif(dim(cancBA)[1],delay,delay+period)
  #   cancBA$yrdx=cancBA$yrdx+pmin(runif(dim(cancBA)[1],delay,delay+period),cancBA$surv)
  cancBAr$seconds=NULL
  cancAr$seconds=NULL
  cancA[cancA$casenum%in%cancAr$casenum,]=cancAr
  #   hist(rbeta(1e4,1,1))
  #    hist(rbeta(1e4,4,4))
  
  
  ######now background B after A cases
  py=cancA$surv
  cancA$seconds=rpois(dim(cancA)[1],(Ab/kb)*(exp(kb*(cancA$age+py))-exp(kb*cancA$age) ) )
  sum(cancA$seconds)
  #   table(cancA$seconds,cut(py,breaks=seq(0,4,0.1)))
  cancA$seqnum[cancA$seconds>0]=1
  cancBA=cancA[cancA$seconds>0,]
  py=py[cancA$seconds>0]
  cancBA$seqnum=2
  cancBA$cancer="B"
  Y=runif(dim(cancBA)[1],Ab*exp(kb*cancBA$age),Ab*exp(kb*(cancBA$age+py)))
  ages=(log(Y)-log(Ab))/kb
  cancBA$yrdx=cancBA$yrdx+ages-cancBA$age
  #   range(ages)
  #   head(cancBA)
  #   cancBA$yrdx=cancBA$yrdx+runif(dim(cancBA)[1],0,py)
  cancBA$seconds=NULL
  cancA$seconds=NULL
  head(cancBA)
  # head(cancA)
  
  
  ######now background B
  B=cbind(popsa[,1:2],cancers=rpois(dim(popsa)[1],Ab*exp(kb*popsa$age)*popsa$py))
  cancB=B[rep(seq_len(nrow(B)), times=B$cancers),]%>%select(-cancers)
  cancB$surv=rexp(dim(cancB)[1],rate=1/tauB)
  cancB$yrdx=cancB$year+runif(dim(cancB)[1],max=0.9999)
  cancB=trimSurv(cancB)
  cancB$cancer="B"
  cancB$trt=sample(trts,dim(cancB)[1],replace=T)
  cancB$casenum=(1e7+1):(1e7+dim(cancB)[1])
  cancB$seqnum=0
  head(cancB)
  
  ######  merge in B
  
  
  canc=rbind(cancA,cancBAr)
  canc=rbind(canc,cancBA)
  canc=rbind(canc,cancB)
  head(canc)
  tail(canc)
  table(canc$seqnum)
  
  
  canc$trt=factor(canc$trt)
  canc$cancer=factor(canc$cancer)
  head(canc)
  cancerS=levels(canc$canc)
  #   sapply(canc,class)
  canc=tbl_df(canc)
  popsa=tbl_df(popsa)
  # and package it all up
  seerSet=list(canc=canc,popsa=popsa,ageStart=min(popsa$age),ageEnd=max(popsa$age),sex="neut",race="neut",
               cancerS=cancerS,yearEnd=max(popsa$year))
  class(seerSet)="seerSet"
  seerSet
} # return a list that can be attached or with-ed in other functions
