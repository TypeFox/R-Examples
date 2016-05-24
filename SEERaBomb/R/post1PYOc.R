post1PYOc=function(canc,brkst=c(0),binIndxt=1,brksy=c(1973),binIndxy=1,brksa=c(0),binIndxa=1,
                   Trt="rad",PYLong=FALSE,yearEnd,firstS,secondS) { 
  # to get rid of check notes. Using list inputs and with will shield these
  surv=yrdx=age=casenum=cancer=trt=yrdx1=seqnum=yrdiff=mn=NULL 
  yrdx2=yrdiffn=cancer1=cancer2=py=year=ageL=ageR=age1=age2=ageM=sem=yrC=ageC=NULL 
  # brkst=c(0);brksy=c(1973,2000);brksa=c(0,50);canc=pf$canc;binIndxt=1;binIndxy=1;binIndxa=1
  if(sum(canc$trt==Trt,na.rm=T)==0) stop(paste0("canc must have a trt column containing",Trt))
  binSt<-levels(cut(brkst+0.1,breaks=c(brkst,100))) #this is just to make a vector of tsd interval/row names 
  binSy<-levels(cut(brksy,breaks=c(brksy,yearEnd+1),right=FALSE,dig.lab=4)) # year at diagnosis break points
  binSa<-levels(cut(brksa+0.1,breaks=c(brksa,126))) # age at diagnosis breaks
  bint=binSt[binIndxt]
  biny=binSy[binIndxy]
  bina=binSa[binIndxa]
  LL=getBinInfo(bint,binSt)["LL"]
  cat("tsd    Int:",bint,"\n")
  cat("yearDx Int:",biny,"\n")
  cat("ageDx  Int:",bina,"\n\n")
  D2=canc%>%filter(seqnum==2) # D2 holds second primaries
  canc$yrC=cut(canc$yrdx,breaks=c(brksy,yearEnd+1),right=FALSE,dig.lab=4)
  canc$ageC=cut(canc$age,breaks=c(brksa,126))
  canc=canc%>%filter(yrC==biny,ageC==bina)
  D0=canc%>%filter(seqnum==0,surv<200,surv>LL,trt==Trt,cancer%in%firstS)
  D0$cancer=factor(D0$cancer,levels=firstS) # get rid of levels not in firstS. 
  # need levels above since apparently pituitary is never irradiated
  D1=canc%>%filter(seqnum==1,trt==Trt,cancer%in%firstS)%>%select(casenum,cancer,yrdx,age,trt)  
  D1$cancer=factor(D1$cancer,levels=firstS) # get rid of levels not in firstS
  D1=D1%>%filter(casenum%in%D2$casenum) # reduce firsts to just those with a second in D2 
  names(D1)[2:5]=c("cancer1","yrdx1","age1","trt1") #rename D1 cols so as not to join by them.
  D2=D2%>%select(casenum,cancer2=cancer,yrdx2=yrdx,age2=age) # reduce D2 to cols we want to slap on 
  D12=left_join(D2,D1,by="casenum") #Keeps all D2 rows, inserts missing where D1 doesn't match. 
  D12=D12%>%filter(!is.na(yrdx1)) # removes firsts before 1973
  D12=D12%>%mutate(yrdiffn=yrdx2-yrdx1)
  D12$yrdiffn[D12$yrdiffn==0]=0.33/12 # if first and second are in the same month, assume 1/3 of a month apart
  D12py=D12 #next 2-lines get cancers observed right, but PY of earlier intervals need to be delivered too
  D12=D12%>%mutate(yrdiff=cut(yrdiffn,breaks=c(-1,brkst,100),include.lowest = TRUE)) 
  D12=D12%>%filter(yrdiff==bint) 
  # I believe multi-arg calls to mutate as below caused R 3.2.0 and 3.2.1 to crash sporadically.
  # PY0=D0%>%mutate(py=getPY(surv,bin,binS,brks),ageL=age+brks[binIndx],year=floor(yrdx+brks[binIndx])) 
  PY0=D0%>%mutate(py=getPY(surv,bint,binSt,brkst)) # getpy leaves zeros when surv end is left of LL
  PY0=PY0%>%filter(py>0)  #get rid of such rows upfront
  PY0=PY0%>%mutate(ageL=age+LL) 
  PY0$year=floor(PY0$yrdx+LL)
  PY0$cancer2="none"
  PY0=PY0%>%select(cancer1=cancer,cancer2,py,ageL,year)
  PY12=D12py%>%mutate(py=getPY(yrdiffn,bint,binSt,brkst)) #getpy leaves zeros when 2nd occurs left of LL
  PY12=PY12%>%filter(py>0)  #get rid of py=0 rows upfront
  PY12=PY12%>%mutate(ageL=age1+LL)
  PY12$year=floor(PY12$yrdx1+LL)
  PY12=PY12%>%select(cancer1,cancer2,py,ageL,year)
  PYL=rbind(PY0,PY12)
  N=dim(PYL)[1]
  binMidPnt=LL+sum(PYL$py)/N/2
  PYL=PYL%>%mutate(ageM=ageL+py/2) 
  PYA=PYL%>%group_by(cancer1)%>%summarize(n=n(),py=sum(py),age=mean(ageM),sig2=var(ageM))%>%mutate(t=LL + py/n/2)
  AgeO=D12%>%filter(cancer2%in%secondS)%>%group_by(cancer1,cancer2)%>%
    summarize(age=mean(age2),sig2=var(age2),n=n()) 
  PYin=PYL%>%select(-cancer1,-cancer2,-ageM)
  LPYin=split(PYin,PYL$cancer1) #splitting done on first cancers, so resulting list names are of first cancers
  LPYinM=lapply(LPYin,as.matrix)
  LPYM=NULL
  # creat a matrix of zeros that is repeatedly the starting point of age-year PY accrual
  yrs=1973:yearEnd
  nyears=length(yrs)
  ages=0.5:125.5
  Zs=matrix(0,ncol=nyears,nrow=length(ages))
  colnames(Zs)=yrs
  rownames(Zs)=ages
  for (i in firstS) {
    PYMat=Zs+0   # fake out system to allocate fresh memory for each instance of the matrix, i.e. each first cancer
    LPYM[[i]]=fillPYM(LPYinM[[i]],PYMat)
  } 
  LD12=split(D12,D12$cancer1) # for getting observed cases in this interval later. Split on first => list names of firsts
  O=t(sapply(LD12,function(x) table(x$cancer2)))
  rownames(O)=names(LD12)
  colnames(O)=levels(D12$cancer2)
  L1=list(LPYM=LPYM,O=O,binMidPnt=binMidPnt,AgeO=AgeO,PYA=PYA) #,PYT=PYT)
  if (PYLong) L1$PYL=PYL # passing the big guy did not slow things down too much
  L1
}



