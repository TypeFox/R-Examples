# simulate.R  (validates tsd implementation using simulated data, makes Figure S4)
rm(list=ls()) 
library(SEERaBomb) 
library(dplyr)  
library(reshape2)
graphics.off()
library(ggplot2)
theme_update(legend.position = c(.75, .815),
             axis.text=element_text(size=rel(1.2)),
             axis.text=element_text(size=rel(1.2)),
             axis.title=element_text(size=rel(1.3)),
             legend.title=element_text(size=rel(1.2)),
             legend.text=element_text(size=rel(1.2)))

simAB<-function(N=1e6,yearEnd=2012,rateA=0.0001,rateB=0.0001,rrAB=5,t1=1,t2=10,mkD=TRUE,te=0) {
  # N=1e6;yearEnd=2012;rateA=0.0001;rateB=0.0001;rrAB=0;library(dplyr)
  lastCN=1 # last casenum
  lastRN=1 # last row number
  brthyrs=1873:2012
  M=matrix(0,nrow=2.3e8*length(brthyrs)*(rateA+rateB),ncol=5)
  colnames(M)<-c("casenum","age","seqnum","yrdx","cancer")
  TC=cbind(matrix(1,ncol=2,nrow=N),matrix(2,ncol=3,nrow=N)) #type of cancer. Need to distinguish 1st from 2nds
  A=rexp(N,rateA)
  B=rexp(N,rateB)
  
  # early detection of background B cases
  BB=B+rexp(N,rateB) # store away now before B gets altered by early detection due to A.
  BE=rexp(N,10) # B Early: on average, latent B are detected within 0.1 yrs of A  
  B=ifelse((B-A<te)&(B-A>0),A+BE,B)
  
  # induction/association of B with A
  if(rrAB>0) {AB=rexp(N,rrAB*rateB) 
              AB=ifelse(AB<t1|AB>t2,Inf,AB) }
  else AB=Inf
  
  AC=cbind(A,A+rexp(N,rateA),A+AB,B,BB)
  system.time(Ord<-t(apply(AC,1,order))) # 7 secs!
  for (i in 1:N) TC[i,]=TC[i,Ord[i,]]
  for (i in 1:N) AC[i,]=AC[i,Ord[i,]]
  I=AC[,1]<100
  TC=TC[I,1:2]
  AC=AC[I,1:2]
  for (k in brthyrs) {
    Age=AC[,1]
    n=length(Age)
    casenumsk=casenums=lastCN:(lastCN+n-1)
    lastCN=lastCN+n
    rownums=lastRN:(lastRN+n-1)
    M[rownums,]=cbind(casenum=casenumsk,age=Age,seqnum=0,yrdx=k+Age,cancer=TC[,1])
    I=AC[,2]<100 
    casenums=casenums[I]
    M[lastRN-1+which(casenumsk%in%casenums),"seqnum"]=1
    lastRN=lastRN+n
    Age=AC[I,2]
    n=length(Age)
    rownums=lastRN:(lastRN+n-1)
    M[rownums,]=cbind(casenum=casenums,age=Age,seqnum=2,yrdx=k+Age,cancer=TC[I,2])
    lastRN=lastRN+n
  } # k loop on birth years
  sum(M[,"casenum"]>0) # 9916+9883 +9948 + 39+ 51 +59
  head(M)
  tail(M)
  dim(M)
  M=M[M[,"casenum"]>0,]
  canc=tbl_df(as.data.frame(M))
  head(canc)
  fixem=canc%>%filter(seqnum==2,yrdx>=(yearEnd+1))
  fixem=fixem$casenum
  canc[canc$casenum%in%fixem,"seqnum"]=0
  canc=canc%>%filter(yrdx>=1973,yrdx<(yearEnd+1))
  canc$surv=100-canc$age
  trimSurv=function(x) {
    x$surv=ifelse(x$surv+x$yrdx>yearEnd+0.9999,yearEnd+0.9999-x$yrdx,x$surv)
    x  }
  canc=trimSurv(canc)
  canc=canc%>%mutate(trt="noRad",year=floor(yrdx))
  canc$trt=factor(canc$trt)
  canc$cancer=factor(canc$cancer,labels=c("A","B"))
  popsa=merge(data.frame(age=0.5:99.5),data.frame(year=1973:yearEnd))
  popsa$py=N
  popsa=tbl_df(popsa)
  cancerS=levels(canc$canc)
  if(mkD) {
    D=canc%>%mutate(age=floor(age)+0.5)
    D=D%>%group_by(cancer,year,age)%>%summarize(cases=n())
    D=cbind(D[,1:3],py=N,D[,4]) 
    D=D%>%mutate(incid=1e5*cases/py,Ecases=ifelse(cancer=="A",rateA,rateB)*1e6,Eincid=Ecases/(N/1e5))
    D=tbl_df(D)
    seerSet=list(canc=canc,popsa=popsa,ageStart=0,ageEnd=100,sex="neut",race="neut",
                 cancerS=cancerS,yearEnd=yearEnd,secondS=c("A","B"),D=D,bfn="nNs0.5e99.5")
  } else {  # if using mk2D downstream
    seerSet=list(canc=canc,popsa=popsa,ageStart=min(popsa$age),ageEnd=max(popsa$age),sex="neut",race="neut",
                 cancerS=cancerS,yearEnd=max(popsa$year))  
  }
  class(seerSet)="seerSet"
  seerSet
} # return a list that can be attached or with-ed in other functions



LO=NULL
LE=NULL
LR=NULL
mybrks=seq(0,20,by=5)
mybrks=c(0,0.75,0.9,1.1,1.25,2,2.5,3,3.5,4,4.75,4.9,5.1,5.25,6)
mybrks=c(0,0.1,0.2,0.3,0.5,0.75,0.9,1.1,1.25,2,2.5,3,3.5,4,4.75,4.9,5.1,5.25,6)
mybrks=c(0,0.1,0.2,0.3,0.5,0.75,1,1.25,2,2.75,3,3.25,4,4.75,5,5.25,6)
for (i in 1:20) {  # i=1 
  print(i)
  n=simAB(rateA=0.001,rateB=0.0002,rrAB=4,t1=1,t2=5,mkD=F,te=0.5)# 100x: 1.011872 1.011866 1.011468 1.012264
  n=mk2D(n)
  n=tsd(n,brks=mybrks,trts=c("noRad")) 
  lab=paste0("b",paste(mybrks,collapse="_"))
  LM=n$L[[lab]]$'noRad'
  LM$Exp
  LM$Obs
  (O=sapply(LM$Obs,cbind)) 
  D=data.frame(first=rep(c("A","B"),2),second=rep(c("A","B"),each=2))
  dO=cbind(D,O)
  (E=sapply(LM$Exp,cbind))
  dE=cbind(D,E)
  LR[[i]]=O/E
  LE[[i]]=E
  LO[[i]]=O
  dO=melt(dO)
  dE=melt(dE)
  M=data.frame(dO,O=dO[,"value"],E=dE[,"value"],t=rep(LM$mids,each=4))
  M=M%>%mutate(RR=O/E)
  g=qplot(x=t,y=RR,col=first,shape=second,data=M,geom=c("line","point"),
          xlab="Years Since First Cancer Diagnosis",ylab="Relative Risk")
  g=g+geom_abline(intercept=1, slope=0)
  print(g)
}

TMP=LR[[1]]
n=length(LR)
for (i in 2:n) TMP=TMP+LR[[i]]
(R=TMP/n)
(R=cbind(D,R))
(LL=qchisq(.025,2*LO[[1]])/(2*LE[[1]]))
(LU=qchisq(.975,2*LO[[1]]+2)/(2*LE[[1]]))
for (i in 2:n) {
  (LL=LL+qchisq(.025,2*LO[[i]])/(2*LE[[i]]))
  (LU=LU+qchisq(.975,2*LO[[i]]+2)/(2*LE[[i]]))
}
(LU=cbind(D,LU=LU/n) )
(LL=cbind(D,LL=LL/n))  
R
dR=melt(R)
dLL=melt(LL)
dLU=melt(LU)

M=data.frame(dR,LL=dLL[,"value"],UL=dLU[,"value"],t=rep(LM$mids,each=4))
head(M)
levels(M$first)
M$first=factor(M$first,labels=c("First Cancer A","First Cancer B"))
theme_update(legend.position = c(.8, .7),
             strip.text=element_text(size=rel(1.2)))
quartz(width=6,height=3.25)
g=qplot(x=t,y=value,col=second,data=M,geom=c("line","point"),xlim=c(0,6),
        xlab="Years Since First Cancer Diagnosis",ylab="Relative Risk")
g=g+facet_grid(~first,scales="free")+geom_abline(intercept=1, slope=0)
g1 <- guide_legend("Second\nCancer")
g=g + guides(color=g1) 
g+  geom_errorbar(aes(ymin=LL,ymax=UL,width=.05))
ggsave("~/Results/amlMDS/ABdipNpeak.png")  
