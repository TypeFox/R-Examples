setClass("uvsd",representation(
                               muN="numeric",
			       alphaN="numeric",
                               betaN="numeric",
                               s2alphaN="numeric",
                               s2betaN="numeric",
                               thetaN="numeric",
                               muS="numeric",
			       alphaS="numeric",
                               betaS="numeric",
                               s2alphaS="numeric",
                               s2betaS="numeric",
                               thetaS="numeric",
                               estN="numeric",
                               estS="numeric",
                               estS2="numeric",
                               estCrit="matrix",
                               blockN="matrix",
                               blockS="matrix",
                               blockS2="matrix",
                               s.crit="array",
                               pD="numeric",
                               DIC="numeric",
                               M="numeric",
                               keep="numeric",
                               b0="matrix",
                               b0S2="numeric",
                               b0Crit="numeric"))

setClass("uvsdSim",representation(
                                  Scond="numeric",
				  cond="numeric",
                                  subj="numeric",
                                  item="numeric",
                                  lag="numeric",
                                  resp="numeric",
                                  muN="numeric",
                                  muS="numeric",
                                  muS2="numeric",
                                  alphaN="numeric",
                                  betaN="numeric",
                                  alphaS="numeric",
                                  betaS="numeric",
                                  alphaS2="numeric",
                                  betaS2="numeric",
				  crit="matrix"))
    


uvsdProbs=function(mean,sd,bounds)
{
cumulative=c(0,pnorm(bounds,mean,sd),1)
p.ij=diff(cumulative)
return(p.ij)
}

uvsdSim=function(NN=2,NS=1,I=30,J=200,K=6,muN=c(-.5,-.2),s2aN=.2,s2bN=.2,muS=.5,s2aS=.2,s2bS=.2,muS2=log(1),s2aS2=0,s2bS2=0,lagEffect=-.001,crit=matrix(rep(c(-1.5,-.5,0,.5,1.5),each=I),ncol=(K-1)))
  {
    #sanity checks:
    if(((J/2)/NN)%%1 !=0 | ((J/2)/NS)%%1 !=0) cat("Number of items must be divisible by number of conditions!\n")
 
    if(length(muN) !=NN) cat("Length of vector of new item means muN must match NN!\n")
    if(length(muS) !=NS) cat("Length of vector of studied item means muS must match NS!\n")

    R=I*J
    alphaN=rnorm(I,0,sqrt(s2aN))
    betaN=rnorm(J,0,sqrt(s2bN))
    alphaS=rnorm(I,0,sqrt(s2aS))
    betaS=rnorm(J,0,sqrt(s2bS))
    alphaS2=rnorm(I,0,sqrt(s2aS2))
    betaS2=rnorm(J,0,sqrt(s2bS2))

#make design matrix 
Scond=cond=subj=item=lag=NULL

for(i in 0:(I-1))
{
tmpScond=sample(rep(0:1,J/2),J)

rep(0:(NN-1),(J/2)/NN)
tmpCond=rep(-99,J)
tmpCond[tmpScond==0]=sample(rep(0:(NN-1),(J/2)/NN),J/2)
tmpCond[tmpScond==1]=sample(rep(0:(NS-1),(J/2)/NS),J/2)

tmpSubj=rep(i,J)
tmpItem=sample(0:(J-1),J)

tmpStudy=sample(tmpItem[tmpScond==1],J/2)
tmpLag=rep(-99,J)
for(j in 1:J) 
{
for(k in 1:(J/2)) 
{
if(tmpStudy[k]==tmpItem[j])
{
tmpLag[j]=((J/2)-k) + j
}}}

Scond=c(Scond,tmpScond)
cond=c(cond,tmpCond)
subj=c(subj,tmpSubj)
item=c(item,tmpItem)
lag=c(lag,tmpLag)
}

lag[Scond==0]=0
lag[Scond==1]=lag[Scond==1]-mean(lag[Scond==1])

    dat=1:R
    for(r in 1:R)
      {
        if(Scond[r]==0) p=uvsdProbs(muN[cond[r]+1]+alphaN[subj[r]+1]+betaN[item[r]+1],1,crit[subj[r]+1,])
        if(Scond[r]==1) p=uvsdProbs(muS[cond[r]+1]+alphaS[subj[r]+1]+betaS[item[r]+1]+lag[r]*lagEffect,sqrt(exp(muS2[cond[r]+1]+alphaS2[subj[r]+1]+betaS2[item[r]+1]+lag[r]*lagEffect)),crit[subj[r]+1,])
        dat[r]=which.max(rmultinom(1,1,p))-1
      }
    
    ret=new("uvsdSim")
    ret@Scond=Scond
    ret@cond=cond
    ret@subj=subj
    ret@item=item
    ret@lag=lag
    ret@resp=dat
    ret@muN=muN
    ret@muS=muS
    ret@muS2=muS2
    ret@alphaN=alphaN
    ret@alphaS=alphaS
    ret@alphaS2=alphaS2
    ret@betaN=betaN
    ret@betaS=betaS
    ret@betaS2=betaS2
    return(ret)
  }


uvsdLogLike=function(R,NN,NS,I,J,K,dat,cond,Scond,subj,item,lag,blockN,blockS,blockS2,crit)
  {
    like=1:R*0
    tmp=.C("logLikeUvsd",as.double(like),as.integer(R),as.integer(NN),as.integer(NS),as.integer(I),as.integer(J),as.integer(dat),as.integer(subj),as.integer(item),as.double(lag),as.integer(cond),as.integer(Scond),as.double(blockN),as.double(blockS),as.double(blockS2),as.double(as.vector(t(crit))),NAOK=TRUE,PACKAGE=.packageName)[[1]]
    return(tmp)
  }


uvsdSample=function(dat,M=10000,keep=(M/10):M,getDIC=TRUE,freeCrit=TRUE,equalVar=FALSE,freeSig2=FALSE,Hier=TRUE,jump=.0001)
{
Scond=dat$Scond
cond=dat$cond
sub=dat$sub
item=dat$item
lag=dat$lag
resp=dat$resp

#DEFINE CONSTANTS
NN=length(levels(as.factor(cond[Scond==0])))
NS=length(levels(as.factor(cond[Scond==1])))
I=length(levels(as.factor(sub)))
J=length(levels(as.factor(item)))
K=length(levels(as.factor(resp)))
ncondN=table(cond[Scond==0])
nsubN=table(sub[Scond==0])
nitemN=table(item[Scond==0])
ncondS=table(cond[Scond==1])
nsubS=table(sub[Scond==1])
nitemS=table(item[Scond==1])
R=dim(dat)[1]
RS=sum(Scond)
RN=R-sum(Scond)
BN=NN+I+J+3
BS=NS+I+J+3

#INDEXING
muN=1:NN
alphaN=(NN+1):(NN+I)
betaN=(NN+I+1):(NN+I+J)
s2alphaN=NN+I+J+1
s2betaN=s2alphaN+1
thetaN=s2betaN+1

muS=1:NS
alphaS=(NS+1):(NS+I)
betaS=(NS+I+1):(NS+I+J)
s2alphaS=NS+I+J+1
s2betaS=s2alphaS+1
thetaS=s2betaS+1

#SPACE AND STARTING VALUES
blockN=matrix(0,nrow=M,ncol=BN)
blockS=blockS2=matrix(0,nrow=M,ncol=BS)
s.crit=array(dim=c(I,7,M))
#blockN[,c(s2alphaN,s2betaN)]=.1
#blockS[,c(s2alphaS,s2betaS)]=.1

if(Hier==0)
  {
    blockN[,c(s2alphaN,s2betaN)]=2
    blockS[,c(s2alphaS,s2betaS)]=2
    blockS2[,c(s2alphaS,s2betaS)]=2
    cat("Non-Hierarchical Model!\n")
    cat("Using Prior Variances of 2.0\n\n")
  }

s.crit[,1,]=-Inf
s.crit[,4,]=0
s.crit[,7,]=Inf
s.crit[,,1]=matrix(rep(c(-Inf,-1.5,-.5,0,.5,1.5,Inf),each=I),ncol=7)
b0=matrix(0,2,2)
met=matrix(c(.01,.01,.001,.001),2,2)
met.crit=rep(.05,I)
b0.crit=rep(0,I)
PropCrit=as.matrix(s.crit[,,1])
if(I==1) PropCrit=t(as.matrix(s.crit[,,1]))
metS2=rep(.3,BS)
b0S2=rep(0,BS)
metS2[thetaS]=.0001

pb=txtProgressBar(min=1,max=M,style=3,width=10)
cat("Starting MCMC\n")
#MCMC LOOP
for(m in 2:M){
#Sample latent data
  wN=rtnorm(RN,getPred(blockN[m-1,],cond[Scond==0],sub[Scond==0],item[Scond==0],lag[Scond==0],NN,I,J,RN),rep(1,RN),s.crit[cbind(sub[Scond==0]+1,resp[Scond==0]+1,m-1)],s.crit[cbind(sub[Scond==0]+1,resp[Scond==0]+2,m-1)])
  wS=rtnorm(RS,getPred(blockS[m-1,],cond[Scond==1],sub[Scond==1],item[Scond==1],lag[Scond==1],NS,I,J,RS),sqrt(exp(getPred(blockS2[m-1,],cond[Scond==1],sub[Scond==1],item[Scond==1],lag[Scond==1],NS,I,J,RS))),s.crit[cbind(sub[Scond==1]+1,resp[Scond==1]+1,m-1)],s.crit[cbind(sub[Scond==1]+1,resp[Scond==1]+2,m-1)])
  
#Sample Blocks
  tmp=sampleNorm(blockN[m-1,],wN,cond[Scond==0],sub[Scond==0],item[Scond==0],lag[Scond==0],NN,I,J,RN,ncondN,nsubN,nitemN,100,.01,.01,met[1,1],met[1,2],1,1,Hier)
  blockN[m,]=tmp[[1]]
  b0[1,]=b0[1,]+tmp[[2]]
  
  if(!freeSig2)
    {
      tmp=sampleNorm(blockS[m-1,],wS,cond[Scond==1],sub[Scond==1],item[Scond==1],lag[Scond==1],NS,I,J,RS,ncondS,nsubS,nitemS,100,.01,.01,met[2,1],met[2,2],exp(blockS2[m-1,muS]),1,Hier)
      blockS[m,]=tmp[[1]]
      b0[2,]=b0[2,]+tmp[[2]]
      
      if(!equalVar)
        {
          blockS2[m,muS]=log(sampleSig2(exp(blockS2[m-1,muS]),blockS[m,],wS,cond[Scond==1],sub[Scond==1],item[Scond==1],lag[Scond==1],NS,ncondS,I,J,.01,.01))
        } 
    }
  
  if(freeSig2)
    {
      tmp=sampleNormb(blockS[m-1,],wS,cond[Scond==1],sub[Scond==1],item[Scond==1],lag[Scond==1],NS,I,J,RS,ncondS,nsubS,nitemS,100,.01,.01,met[2,1],met[2,2],blockS2[m-1,],1,Hier)
      blockS[m,]=tmp[[1]]
      b0[2,]=b0[2,]+tmp[[2]]
      
      tmp=sampleSig2b(blockS2[m-1,],wS,cond[Scond==1],sub[Scond==1],item[Scond==1],lag[Scond==1],NS,I,J,RS,ncondS,nsubS,nitemS,100,.01,.01,metS2,blockS[m,],1,Hier)
      blockS2[m,]=tmp[[1]]
      b0S2=b0S2+tmp[[2]]
    }
   
#Sample Criteria
  if(!freeCrit)
    {
      PropCrit[,c(2,3,5,6)]=s.crit[,c(2,3,5,6),m-1]+rep(rnorm(4,0,met.crit[1]),each=I)
      violate=pmin(1,sum(diff(PropCrit[1,2:6])<0))
      if(violate) PropCrit=t(as.matrix(s.crit[,,m-1]))
      likeOld=sum(uvsdLogLike(R,NN,NS,I,J,K,resp,cond,Scond,sub,item,lag,blockN[m,],blockS[m,],blockS2[m,],s.crit[,,m-1]))
      likeProp=sum(uvsdLogLike(R,NN,NS,I,J,K,resp,cond,Scond,sub,item,lag,blockN[m,],blockS[m,],blockS2[m,],PropCrit))
      accept=rbinom(1,1,min(1,(1-violate)*exp(likeProp-likeOld)))
      b0.crit=b0.crit+accept
      s.crit[,c(2,3,5,6),m]=accept*PropCrit[,c(2,3,5,6)] + (1-accept)*s.crit[,c(2,3,5,6),m-1]
    }
  
  if(freeCrit)
    {
      PropCrit[,c(2,3,5,6)]=s.crit[,c(2,3,5,6),m-1]+matrix(rnorm(I*4,0,rep(met.crit,4)),ncol=4)
      violate=pmin(1,colSums((apply(PropCrit[,c(2,3,4,5,6)],1,diff))<0))
      PropCrit[violate==1,]=s.crit[violate==1,,m-1] 
      likeOld=tapply(uvsdLogLike(R,NN,NS,I,J,K,resp,cond,Scond,sub,item,lag,blockN[m,],blockS[m,],blockS2[m,],s.crit[,,m-1]),sub,sum)
      likeProp=tapply(uvsdLogLike(R,NN,NS,I,J,K,resp,cond,Scond,sub,item,lag,blockN[m,],blockS[m,],blockS2[m,],PropCrit),sub,sum)
      accept=rbinom(I,1,pmin(1,(1-violate)*exp(likeProp-likeOld)))
      
      b0.crit=b0.crit+accept
      s.crit[accept==1,,m]=PropCrit[accept==1,]
      s.crit[accept==0,,m]=s.crit[accept==0,,m-1]
    }
  
#AUTOTUNING
  if(m>20 & m<keep[1] & m%%5==0)
    {
      met=met+(b0/m<.3)*matrix(-jump,2,2) +(b0/m>.5)*matrix(jump,2,2)
      met.crit=met.crit+((b0.crit/m<.3)*-jump + (b0.crit/m>.5)*jump)
      metS2=metS2+((b0S2/m<.3)*-jump + (b0S2/m>.5)*jump)
      
      met[met<jump]=jump
      met.crit[met.crit<jump]=jump
      metS2[metS2<jump]=jump
    }
  
  if(m==keep[1])
    {
      cat("\n\nContinuing with following acceptance probabilities:\n")
      cat(round(b0/m,2),"\n",fill=TRUE)
      if(freeSig2) cat(round(b0S2/m,2),"\n",fill=TRUE)
      cat(round(b0.crit/m,2),"\n",fill=TRUE)
      cat("If they are <.2 or >.6, adjust jump or keep\n")
      
      #cat("MetS2: \n",metS2,fill=TRUE)
    }
  setTxtProgressBar(pb, m)
}

close(pb)

estN=colMeans(blockN[keep,])
estS=colMeans(blockS[keep,])
estS2=colMeans(blockS2[keep,])
estCrit=apply(s.crit[,,keep],c(1,2),mean)
if(I==1) estCrit=as.matrix(rowMeans(s.crit[1,,keep]))


#GET DIC
DIC=pD=NA
if(getDIC)
{
cat("Getting DIC\n")
pb=txtProgressBar(min=keep[1],max=tail(keep,1),style=3,width=10)

D0=0
for(m in keep)
  {
    D0=D0+sum(-2*uvsdLogLike(R,NN,NS,I,J,K,resp,cond,Scond,sub,item,lag,blockN[m,],blockS[m,],blockS2[m,],s.crit[,,m]))
    setTxtProgressBar(pb, m)
  }
D0=D0/length(keep)
Dhat=sum(-2*uvsdLogLike(R,NN,NS,I,J,K,resp,cond,Scond,sub,item,lag,estN,estS,estS2,estCrit))
pD=D0-Dhat
DIC=pD+D0
}


u=new("uvsd")
u@muN=muN
u@alphaN=alphaN
u@betaN=betaN
u@s2alphaN=s2alphaN
u@s2betaN=s2betaN
u@thetaN=thetaN
u@muS=muS
u@alphaS=alphaS
u@betaS=betaS
u@s2alphaS=s2alphaS
u@s2betaS=s2betaS
u@thetaS=thetaS
u@estN=estN
u@estS=estS
u@estS2=estS2
u@estCrit=estCrit
u@blockN=blockN[keep,]
u@blockS=blockS[keep,]
u@blockS2=blockS2[keep,]
u@s.crit=s.crit[,,keep]
u@pD=pD
u@DIC=DIC
u@M=M
u@keep=keep
u@b0=b0/M
u@b0S2=b0S2/M
u@b0Crit=b0.crit/M

cat("Finished \n\n")
return(u)
}
