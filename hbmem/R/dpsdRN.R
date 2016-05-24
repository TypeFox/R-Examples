dpsdRNSim=function(NN=2,NS=1,I=30,J=200,K=6,muN=c(-.7,-.5),s2aN=.2,s2bN=.2,muS=0,s2aS=.2,s2bS=.2,muR=qnorm(.25),phiA=-1,phiB=-1,crit=matrix(rep(c(-1.6,-.5,0,.5,1.6),each=I),ncol=(K-1)))
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
    alphaR=phiA*alphaN
    betaR=phiB*betaN
 
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

    dat=1:R
    for(r in 1:R)
      {
        if(Scond[r]==0) p=dpsdProbs(0,muN[cond[r]+1]+alphaN[subj[r]+1]+betaN[item[r]+1],crit[subj[r]+1,])
        if(Scond[r]==1) p=dpsdProbs(pnorm(muR[cond[r]+1]+alphaR[subj[r]+1]+betaR[item[r]+1]),muS[cond[r]+1]+alphaS[subj[r]+1]+betaS[item[r]+1],crit[subj[r]+1,])
        dat[r]=which.max(rmultinom(1,1,p))-1
      }

     
    ret=new("dpsdSim")
    ret@Scond=Scond
    ret@cond=cond
    ret@subj=subj
    ret@item=item
    ret@lag=lag
    ret@resp=dat
    ret@muN=muN
    ret@muS=muS
    ret@muR=muR
    ret@alphaN=alphaN
    ret@alphaS=alphaS
    ret@alphaR=alphaR
    ret@betaN=betaN
    ret@betaS=betaS
    ret@betaR=betaR
    return(ret)
  }


dpsdRNSample=function(dat,M=5000,keep=(M/10):M,getDIC=TRUE,jump=.001)
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
RN=R-RS
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
blockS=blockR=matrix(0,nrow=M,ncol=BS)
blockN[1,muN]=rep(-.2,NN)
blockN[1,c(s2alphaN,s2betaN)]=.1
blockS[1,c(s2alphaS,s2betaS)]=.1
phiA=phiB=1:M*0
phiA[1]=phiB[1]=-1
blockR[1,c(s2alphaS,s2betaS)]=c(phiA[1],phiB[1])
wS=wR=rnorm(RS)
s.crit=array(dim=c(I,7,M))
s.crit[,1,]=-Inf
s.crit[,4,]=0
s.crit[,7,]=Inf
s.crit[,,1]=matrix(rep(c(-Inf,-1,-.5,0,.5,1,Inf),each=I),ncol=7)
b0=matrix(0,3,2)
met=matrix(rep(.01,6),3,2)
met.crit=rep(.05,I)
b0.crit=rep(0,I)
PropCrit=s.crit[,,1]

notk=Scond==1 & resp<(K-1)
Rnotk=sum(notk)
isk= Scond==1 & resp==(K-1) 
Risk=sum(isk)

pb=txtProgressBar(min=1,max=M,style=3,width=10)
cat("Starting MCMC\n")
#MCMC LOOP
for(m in 2:M){
#Sample latent data
wN=rtnorm(RN,getPred(blockN[m-1,],cond[Scond==0],sub[Scond==0],item[Scond==0],lag[Scond==0],NN,I,J,RN),rep(1,RN),s.crit[cbind(sub[Scond==0]+1,resp[Scond==0]+1,m-1)],s.crit[cbind(sub[Scond==0]+1,resp[Scond==0]+2,m-1)])

wS[notk[Scond==1]]=rtnorm(Rnotk,getPred(blockS[m-1,],cond[notk],sub[notk],item[notk],lag[notk],NS,I,J,Rnotk),rep(1,Rnotk),s.crit[cbind(sub[notk]+1,resp[notk]+1,m-1)],s.crit[cbind(sub[notk]+1,resp[notk]+2,m-1)])
wR[notk[Scond==1]]=rtnorm(Rnotk,getPred(blockR[m-1,],cond[notk],sub[notk],item[notk],lag[notk],NS,I,J,Rnotk),rep(1,Rnotk),rep(-Inf,Rnotk),rep(0,Rnotk))

ind=isk[Scond==1] & wR>0
wS[ind]=rnorm(sum(ind),getPred(blockS[m-1,],cond[Scond==1][ind],sub[Scond==1][ind],item[Scond==1][ind],lag[Scond==1][ind],NS,I,J,sum(ind)),1)
ind=isk[Scond==1] & wR<0
wS[ind]=rtnorm(sum(ind),getPred(blockS[m-1,],cond[Scond==1][ind],sub[Scond==1][ind],item[Scond==1][ind],lag[Scond==1][ind],NS,I,J,sum(ind)),rep(1,sum(ind)),s.crit[cbind(sub[Scond==1][ind]+1,K,m-1)],rep(Inf,sum(ind)))

ind=isk[Scond==1] & wS>s.crit[cbind(sub[Scond==1]+1,K,m-1)]
wR[ind]=rnorm(sum(ind),getPred(blockR[m-1,],cond[Scond==1][ind],sub[Scond==1][ind],item[Scond==1][ind],lag[Scond==1][ind],NS,I,J,sum(ind)),1)
ind=isk[Scond==1] & wS<s.crit[cbind(sub[Scond==1]+1,K,m-1)]
wR[ind]=rtnorm(sum(ind),getPred(blockR[m-1,],cond[Scond==1][ind],sub[Scond==1][ind],item[Scond==1][ind],lag[Scond==1][ind],NS,I,J,sum(ind)),rep(1,sum(ind)),rep(0,sum(ind)),rep(Inf,sum(ind)))

#Sample Studied Block 
tmp=sampleNorm(blockS[m-1,],wS,cond[Scond==1],sub[Scond==1],item[Scond==1],lag[Scond==1],NS,I,J,RS,ncondS,nsubS,nitemS,100,.01,.01,met[2,1],met[2,2],1,1)
blockS[m,]=tmp[[1]]
b0[2,]=b0[2,]+tmp[[2]]

#Sample muR
tmp=wR-blockR[m-1,alphaS[sub[Scond==1]+1]]-blockR[m-1,betaS[item[Scond==1]+1]]
tmp=tapply(tmp,cond[Scond==1],sum)
a=1/(ncondS + 1/100);
blockR[m,muS]=rnorm(NS,tmp*a,sqrt(a))

#Sample muN
tmp=wN-blockN[m-1,alphaN[sub[Scond==0]+1]]-blockN[m-1,betaN[item[Scond==0]+1]]
tmp=tapply(tmp,cond[Scond==0],sum)
a=1/(ncondN + 1/100);
blockN[m,muN]=rnorm(NN,tmp*a,sqrt(a))

#Sample AlphaN
tmpN=wN-blockN[m,muN[cond[Scond==0]+1]]-blockN[m-1,betaN[item[Scond==0]+1]]
tmpN=tapply(tmpN,sub[Scond==0],sum)
tmpR=phiA[m-1]*(wR-blockR[m,muS[cond[Scond==1]+1]]-blockR[m-1,betaS[item[Scond==1]+1]])
tmpR=tapply(tmpR,sub[Scond==1],sum)
tmp=tmpN+tmpR
a=1/(nsubN + nsubS*phiA[m-1]^2+ 1/blockN[m-1,s2alphaN])
blockN[m,alphaN]=rnorm(I,tmp*a,sqrt(a))

#Sample s2AlphaN
postB=sum(blockN[m,alphaN]^2)/2 + .01
blockN[m,s2alphaN]=1/rgamma(1,I/2 + .01,rate=postB)

#Sample BetaN
tmpN=wN-blockN[m,muN[cond[Scond==0]+1]]-blockN[m-1,alphaN[sub[Scond==0]+1]]
tmpN=tapply(tmpN,item[Scond==0],sum)
tmpR=phiB[m-1]*(wR-blockR[m,muS[cond[Scond==1]+1]]-blockR[m-1,alphaS[sub[Scond==1]+1]])
tmpR=tapply(tmpR,item[Scond==1],sum)
tmp=tmpN+tmpR
a=1/(nitemN + nitemS*phiB[m-1]^2 + 1/blockN[m-1,s2betaN])
blockN[m,betaN]=rnorm(J,tmp*a,sqrt(a))

#Sample S2BetaN
postB=sum(blockN[m,betaN]^2)/2 + .01
blockN[m,s2betaN]=1/rgamma(1,J/2 + .01,rate=postB)

#Sample PhiAlpha
tmp=sum(blockN[m-1,alphaN[sub[Scond==1]+1]]*
(wR-blockR[m,muS[cond[Scond==1]+1]]-phiB[m-1]*blockN[m-1,betaN[item[Scond==1]+1]]))
a=1/(sum(nsubS*blockN[m-1,alphaN]^2) + 1/100);
phiA[m]=rnorm(1,tmp*a,sqrt(a))

#Sample PhiBeta
tmp=sum(blockN[m-1,betaN[item[Scond==1]+1]]*
(wR-blockR[m,muS[cond[Scond==1]+1]]-phiA[m]*blockN[m-1,alphaN[sub[Scond==1]+1]]))
a=1/(sum(nitemS*blockN[m-1,betaN]^2) + 1/100);
phiB[m]=rnorm(1,tmp*a,sqrt(a))

#Set BlockR to BlockN for convience
blockR[m,alphaS]=phiA[m]*blockN[m,alphaN]
blockR[m,betaS]=phiB[m]*blockN[m,betaN]
blockR[m,thetaS]=0
blockR[m,s2alphaS]=phiA[m]
blockR[m,s2betaS]=phiB[m]



#Sample Criteria
PropCrit[,c(2,3,5,6)]=s.crit[,c(2,3,5,6),m-1]+rnorm(I*4,0,rep(met.crit,4))
violate=pmin(1,colSums((apply(PropCrit[,c(2,3,4,5,6)],1,diff))<0))
PropCrit[violate==1,]=s.crit[violate==1,,m-1]

likeOld=tapply(dpsdLogLike(R,NN,NS,I,J,K,resp,cond,Scond,sub,item,lag,blockN[m,],blockS[m,],blockR[m,],s.crit[,,m-1]),sub,sum)
likeProp=tapply(dpsdLogLike(R,NN,NS,I,J,K,resp,cond,Scond,sub,item,lag,blockN[m,],blockS[m,],blockR[m,],PropCrit),sub,sum)
accept=rbinom(I,1,pmin(1,(1-violate)*exp(likeProp-likeOld)))
b0.crit=b0.crit+accept
s.crit[accept==1,,m]=PropCrit[accept==1,]
s.crit[accept==0,,m]=s.crit[accept==0,,m-1]


#AUTOTUNING
if(m>20 & m<keep[1] & m%%10==0)
  {
    met=met+(b0/m<.3)*matrix(-jump,3,2) +(b0/m>.5)*matrix(jump,3,2)
    met.crit=met.crit+((b0.crit/m<.3)*-jump + (b0.crit/m>.5)*jump)
    met[met<jump]=jump
    met.crit[met.crit<jump]=jump
  }

if(m==keep[1])
  {
    cat("\n\nContinuing with following acceptance probabilities:\n")
    cat(round(b0/m,2),"\n",fill=TRUE)
    cat(round(b0.crit/m,2),"\n",fill=TRUE)
    cat("If they are <.2 or >.6, adjust jump or keep\n")
  }

setTxtProgressBar(pb, m)
}
close(pb)

estN=colMeans(blockN[keep,])
estS=colMeans(blockS[keep,])
estR=colMeans(blockR[keep,])
estCrit=apply(s.crit[,,keep],c(1,2),mean)

#GET DIC
DIC=pD=NA
if(getDIC)
{
print("Getting DIC")
pb=txtProgressBar(min=keep[1],max=tail(keep,1),style=3,width=10)

D0=0
for(m in keep)
  {
    D0=D0+sum(-2*dpsdLogLike(R,NN,NS,I,J,K,resp,cond,Scond,sub,item,lag,blockN[m,],blockS
      [m,],blockR[m,],s.crit[,,m]))
    setTxtProgressBar(pb, m)
  }
D0=D0/length(keep)
Dhat=sum(-2*dpsdLogLike(R,NN,NS,I,J,K,resp,cond,Scond,sub,item,lag,estN,estS,estR,estCrit))
pD=D0-Dhat
DIC=pD+D0
}


u=new("dpsd")
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
u@estR=estR
u@estCrit=estCrit
u@blockN=blockN[keep,]
u@blockS=blockS[keep,]
u@blockR=blockR[keep,]
u@s.crit=s.crit[,,keep]
u@pD=pD
u@DIC=DIC
u@M=M
u@keep=keep
u@b0=b0/M
u@b0Crit=b0.crit/M
print("")
print("Finished")
return(u)
}
