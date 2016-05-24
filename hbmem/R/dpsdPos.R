dpsdPosSim=function(NN=1,NS=2,I=30,J=200,K=6,muN=-.7,s2aN=.2,s2bN=.2,muD=c(0,.5),s2aD=.2,s2bD=.2,muR=qnorm(c(.2,.4)),s2aR=.2,s2bR=.2,crit=matrix(rep(c(-1.6,-.5,0,.5,1.6),each=I),ncol=(K-1)))
  {
    R=I*J
    alphaN=rnorm(I,0,sqrt(s2aN))
    betaN=rnorm(J,0,sqrt(s2bN))
    alphaD=rnorm(I,0,sqrt(s2aD))
    betaD=rnorm(J,0,sqrt(s2bD))
    alphaR=rnorm(I,0,sqrt(s2aR))
    betaR=rnorm(J,0,sqrt(s2bR))
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
        if(Scond[r]==1) p=dpsdProbs(pnorm(muR[cond[r]+1]+alphaR[subj[r]+1]+betaR[item[r]+1]),muN[1]+alphaN[subj[r]+1]+betaN[item[r]+1]+exp(muD[cond[r]+1]+alphaD[subj[r]+1]+betaD[item[r]+1]),crit[subj[r]+1,])
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
    ret@muS=muD
    ret@muR=muR
    ret@alphaN=alphaN
    ret@alphaS=alphaD
    ret@alphaR=alphaR
    ret@betaN=betaN
    ret@betaS=betaD
    ret@betaR=betaR
    return(ret)
  }

dpsdPosLogLike=function (R, NN, NS, I, JN, JS, K, dat, cond, Scond, sub, item, 
    lag, blockN, blockD, blockR, crit) 
{
    l = 1:R * 0
    .C("logLikeDpsd", as.double(l), as.integer(R), as.integer(NN), 
        as.integer(NS), as.integer(I), as.integer(JN), as.integer(JS), 
        as.integer(K), as.integer(dat), as.integer(cond), as.integer(Scond), 
        as.integer(sub), as.integer(item), as.double(lag), as.double(blockN), 
        as.double(blockD), as.double(blockR), as.double(as.vector(t(crit))), 
        NAOK = TRUE, PACKAGE = .packageName)[[1]]
}

dpsdPosSample=function(dat,M=5000,keep=(M/10):M,getDIC=TRUE,jump=.01)
{
Scond=dat$Scond
cond=dat$cond
sub=dat$sub
item=dat$item
lag=dat$lag
resp=dat$resp

#in case items don't appear in new and studied conditions
tmpN=as.numeric(as.factor(item[Scond==0]))
tmpS=as.numeric(as.factor(item[Scond==1]))
item[Scond==0]=tmpN-1
item[Scond==1]=tmpS-1

#DEFINE CONSTANTS
NN=length(levels(as.factor(cond[Scond==0])))
if(NN!=1)
{
	print("ERROR: Can only have 1 baseline condition")
	break()
}
condN=rep(0,length(Scond))
NS=length(levels(as.factor(cond[Scond==1])))
I=length(levels(as.factor(sub)))
JN=length(levels(as.factor(item[Scond==0])))
JS=length(levels(as.factor(item[Scond==1])))
J=max(c(JN,JS))
K=length(levels(as.factor(resp)))
ncondN=table(cond[Scond==0])
nsubT=table(sub)
nitemT=table(item)
ncondS=table(cond[Scond==1])
nsubS=table(sub[Scond==1])
nitemS=table(item[Scond==1])
R=dim(dat)[1]
RS=sum(Scond)
RN=R-RS
BN=NN+I+JN+3
BS=NS+I+JS+3

#INDEXING
muN=1:NN
alphaN=(NN+1):(NN+I)
betaN=(NN+I+1):(NN+I+JN)
s2alphaN=NN+I+JN+1
s2betaN=s2alphaN+1
thetaN=s2betaN+1
muS=1:NS
alphaS=(NS+1):(NS+I)
betaS=(NS+I+1):(NS+I+JS)
s2alphaS=NS+I+JS+1
s2betaS=s2alphaS+1
thetaS=s2betaS+1

#SPACE AND STARTING VALUES
blockN=matrix(0,nrow=M,ncol=BN)
blockD=blockR=matrix(0,nrow=M,ncol=BS)
wS=wR=rnorm(RS)
allwN=rnorm(R)
s.crit=array(0,dim=c(I,K+1,M))
s.crit[,1,]=-Inf
s.crit[,K+1,]=Inf
s.crit[,,1]=matrix(rep(c(-Inf,-1,-.5,0,.5,1,Inf),each=I),ncol=7)

b0=matrix(0,2,2)
met=matrix(.01,2,2)
b0D=rep(0,BS)
metD=c(.01,rep(.2,(BS-4)),.05,.05,.01)
met.crit=rep(.05,I)
b0.crit=rep(0,I)
PropCrit=s.crit[,,1]

notk=Scond==1 & resp<(K-1)
Rnotk=sum(notk)
isk= Scond==1 & resp==(K-1) 
Risk=sum(isk)

pb=txtProgressBar(min=1,max=M,style=3,width=10)
print("Starting MCMC")
#MCMC LOOP
for(m in 2:M){
#Sample latent data
meanN=getPred(blockN[m-1,],rep(0,length(sub)),sub,item,lag,1,I,J,R)
meanD=exp(getPred(blockD[m-1,],cond[Scond==1],sub[Scond==1],item[Scond==1],lag[Scond==1],NS,I,JS,RS))
meanS=meanN[Scond==1]+meanD
meanR=getPred(blockR[m-1,],cond[Scond==1],sub[Scond==1],item[Scond==1],lag[Scond==1],NS,I,JS,RS)

wN=rtnorm(RN,meanN[Scond==0],rep(1,RN),s.crit[cbind(sub[Scond==0]+1,resp[Scond==0]+1,m-1)],s.crit[cbind(sub[Scond==0]+1,resp[Scond==0]+2,m-1)])

wS[notk[Scond==1]]=rtnorm(Rnotk,meanS[notk],rep(1,Rnotk),s.crit[cbind(sub[notk]+1,resp[notk]+1,m-1)],s.crit[cbind(sub[notk]+1,resp[notk]+2,m-1)])
wR[notk[Scond==1]]=rtnorm(Rnotk,meanR[notk],rep(1,Rnotk),rep(-Inf,Rnotk),rep(0,Rnotk))

ind=isk[Scond==1] & wR>0
wS[ind]=rnorm(sum(ind),meanS[ind],1)
ind=isk[Scond==1] & wR<0
wS[ind]=rtnorm(sum(ind),meanS[ind],rep(1,sum(ind)),s.crit[cbind(sub[ind]+1,K,m-1)],rep(Inf,sum(ind)))

ind=isk[Scond==1] & wS>s.crit[cbind(sub[Scond==1]+1,K,m-1)]
wR[ind]=rnorm(sum(ind),meanR[ind],1)
ind=isk[Scond==1] & wS<s.crit[cbind(sub[Scond==1]+1,K,m-1)]
wR[ind]=rtnorm(sum(ind),meanR[ind],rep(1,sum(ind)),rep(0,sum(ind)),rep(Inf,sum(ind)))

###MIKE LEFT OFF HERE TO GO SEE TRON###

#Sample Blocks
wD=wS-meanN[Scond==1]
tmp=samplePosNorm(blockD[m-1,],wD,sub[cond==1],item[cond==1],lag[cond==1],I,J,RS,10,.01,.01,metD,1,1)
blockD[m,]=tmp[[1]]
b0D=b0D+tmp[[2]]

meanD=exp(getPred(blockD[m,],sub,item,lag,I,J,R))
allwN[cond==0]=wN
allwN[cond==1]=wS-meanD[cond==1]
tmp=sampleNorm(blockN[m-1,],allwN,sub,item,lag,I,J,R,nsubT,nitemT,10,.01,.01,met[1,1],met[1,2],1,0)
blockN[m,]=tmp[[1]]
b0[1,]=b0[1,]+tmp[[2]]

tmp=sampleNorm(blockR[m-1,],wR,sub[cond==1],item[cond==1],lag[cond==1],I,J,RS,nsubS,nitemS,10,.01,.01,met[2,1],met[2,2],1,1)
blockR[m,]=tmp[[1]]
b0[2,]=b0[2,]+tmp[[2]]

#Sample Criteria
PropCrit[,c(2,3,5,6)]=s.crit[,c(2,3,5,6),m-1]+rnorm(I*4,0,rep(met.crit,4))
violate=pmin(1,colSums((apply(PropCrit[,c(2,3,4,5,6)],1,diff))<0))
PropCrit[violate==1,]=s.crit[violate==1,,m-1]

likeOld=tapply(dpsdPosLogLike(R,I,J,K,resp,cond,sub,item,lag,blockN[m,],blockD[m,],blockR[m,],s.crit[,,m-1]),sub,sum)
likeProp=tapply(dpsdPosLogLike(R,I,J,K,resp,cond,sub,item,lag,blockN[m,],blockD[m,],blockR[m,],PropCrit),sub,sum)
accept=rbinom(I,1,pmin(1,(1-violate)*exp(likeProp-likeOld)))
b0.crit=b0.crit+accept
s.crit[accept==1,,m]=PropCrit[accept==1,]
s.crit[accept==0,,m]=s.crit[accept==0,,m-1]


#AUTOTUNING
if(m>20 & m<keep[1] & m%%10==0)
  {
    met=met+(b0/m<.3)*matrix(-jump,2,2) +(b0/m>.5)*matrix(jump,2,2)
    metD=metD+(b0D/m<.3)*-jump + (b0D/m>.5)*jump
    met.crit=met.crit+((b0.crit/m<.3)*-jump + (b0.crit/m>.5)*jump)
    met[met<jump]=jump
    metD[metD<jump]=jump
    met.crit[met.crit<jump]=jump
  }

if(m==keep[1])
  {
    print("")
    print("Continuning with following acceptance probabilities:")
    print(b0/m)
    print(b0D/m)
    print(b0.crit/m)
    print("If they are <.2 or >.6, adjust jump or keep")
  }

setTxtProgressBar(pb, m)
}
close(pb)

estN=colMeans(blockN[keep,])
estD=colMeans(blockD[keep,])
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
    D0=D0+sum(-2*dpsdPosLogLike(R,I,J,K,resp,cond,sub,item,lag,blockN[m,],blockD[m,],blockR[m,],s.crit[,,m]))
    setTxtProgressBar(pb, m)
  }
D0=D0/length(keep)
Dhat=sum(-2*dpsdPosLogLike(R,I,J,K,resp,cond,sub,item,lag,estN,estD,estR,estCrit))
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
u@estS=estD
u@estR=estR
u@estCrit=estCrit
u@blockN=blockN[keep,]
u@blockS=blockD[keep,]
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
