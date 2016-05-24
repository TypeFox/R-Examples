dpsdFullSim=function(I=30,J=200,K=6,muN=-.7,s2aN=.2,s2bN=.2,muD=0,s2aD=.2,s2bD=.2,muR=qnorm(.25),phiA=1,etaA=0,phiB=1,etaB=0,crit=matrix(rep(c(-1.6,-.5,0,.5,1.6),each=I),ncol=(K-1)))
  {
    R=I*J
    alphaN=rnorm(I,0,sqrt(s2aN))
    betaN=rnorm(J,0,sqrt(s2bN))
    alphaD=rnorm(I,0,sqrt(s2aD))
    betaD=rnorm(J,0,sqrt(s2bD))
    alphaR=rnorm(I,phiA*alphaD,sqrt(etaA))
    betaR=rnorm(J,phiB*betaD,sqrt(etaB))
    subj=rep(0:(I-1),each=J)
    item=rep(0:(J-1),I)
    lag=rep(0,R)
    cond.sub.A=rep(0:1,J/2)
    cond.sub.B=rep(1:0,J/2)
    cond=rep(c(cond.sub.A,cond.sub.B),I/2)

    dat=1:R
    for(r in 1:R)
      {
        meanN=muN+alphaN[subj[r]+1]+betaN[item[r]+1]
        if(cond[r]==0) p=dpsdProbs(0,meanN,crit[subj[r]+1,])
        if(cond[r]==1)
          {
            delta=exp(muD + alphaD[subj[r]+1]+betaD[item[r]+1])
            meanS=meanN+delta
            pR=pnorm(muR+alphaR[subj[r]+1]+betaR[item[r]+1])
            p=dpsdProbs(pR,meanS,crit[subj[r]+1,])
          }
           dat[r]=which.max(rmultinom(1,1,p))-1
      }

    ret=new("dpsdSim")
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

cond.muD=function(x,mu,alpha,beta,theta,sub,item,lag,sig2mu)
{
mean=exp(mu+alpha[sub+1]+beta[item+1]+theta*lag)
return(-.5*sum(mean^2 - 2*x*mean) + (mu^2/sig2mu))
}

cond.alphaD=function(wD,muD,alphaD,betaD,thetaD,alphaR,etaA,phiA,sub,item,lag,sig2alpha)
{
meanD=exp(muD+alphaD[sub+1]+betaD[item+1]+thetaD*lag)
like=-.5*(meanD^2 - 2*wD*meanD)
return(tapply(like,sub,sum) - .5*(alphaD^2/sig2alpha) - ((phiA*alphaD)^2 - 2*phiA*alphaD*alphaR)/(2*etaA))
}

cond.betaD=function(wD,muD,alphaD,betaD,thetaD,betaR,etaB,phiB,sub,item,lag,sig2beta)
{
meanD=exp(muD+alphaD[sub+1]+betaD[item+1]+thetaD*lag)
like=-.5*(meanD^2 - 2*wD*meanD) 
return(tapply(like,item,sum) - .5*(betaD^2/sig2beta) - ((phiB*betaD)^2 - 2*phiB*betaD*betaR)/(2*etaB))
}

cond.thetaD=function(wD,muD,alphaD,betaD,thetaD,wR,muR,alphaR,betaR,sub,item,lag,sig2theta)
{
meanD=exp(muD+alphaD[sub+1]+betaD[item+1]+thetaD*lag)
likeD=-.5*(meanD^2 - 2*wD*meanD) 
likeR=-.5*((thetaD*lag)^2 - 2*thetaD*lag*(wR-muR-alphaR[sub+1]-betaR[item+1]))
return(sum(likeD)+sum(likeR)-.5*thetaD^2/sig2theta)
}

dpsdFullSample=function(dat,M=5000,keep=(M/10):M,getDIC=TRUE,jump=.001)
{
cond=dat$cond
sub=dat$sub
item=dat$item
lag=dat$lag
resp=dat$resp


#DEFINE CONSTANTS
I=length(levels(as.factor(sub)))
J=length(levels(as.factor(item)))
K=length(levels(as.factor(resp)))
nsubT=table(sub)
nitemT=table(item)
nsubS=table(sub[cond==1])
nitemS=table(item[cond==1])
R=dim(dat)[1]
RS=sum(cond)
RN=R-sum(cond)
B=I+J+4
sig2mu=100
a=b=.01

#INDEXING
mu=1
alpha=2:(I+1)
beta=(max(alpha)+1):(I+J+1)
s2alpha=max(beta)+1
s2beta=s2alpha+1
theta=s2beta+1

#SPACE AND STARTING VALUES
NK=length(keep)
blockN=blockD=blockR=matrix(0,nrow=2,ncol=B)
phi=rep(0,2)
s.crit=array(dim=c(I,7,2))
save.blockN=save.blockD=save.blockR=matrix(NA,nrow=NK,ncol=B)
save.crit=array(dim=c(I,7,NK))
save.phi=rep(NA,NK)
blockN[1,s2alpha]=blockN[1,s2beta]=blockD[1,s2alpha]=blockD[1,s2beta]=blockR[1,s2alpha]=blockR[1,s2beta]=.001
allwN=rnorm(R)
wS=wR=rnorm(RS)
s.crit[,1,]=-Inf
s.crit[,4,]=0
s.crit[,7,]=Inf
s.crit[,,1]=matrix(rep(c(-Inf,-1,-.5,0,.5,1,Inf),each=I),ncol=7)
b0=matrix(0,2,2)
met=matrix(.01,2,2)
b0D=rep(0,B)
metD=c(.01,rep(.2,(B-4)),.05,.05,.01)
met.crit=rep(.05,I)
b0.crit=rep(0,I)
PropCrit=s.crit[,,1]

notk=cond==1 & resp<(K-1)
Rnotk=sum(notk)
isk= cond==1 & resp==(K-1) 
Risk=sum(isk)

pb=txtProgressBar(min=1,max=M,style=3,width=10)
print("Starting MCMC")
#MCMC LOOP
for(m in 2:M){
#Sample latent data
meanN=getPred(blockN[1,],sub,item,lag,I,J,R)
meanD=exp(getPred(blockD[1,],sub,item,lag,I,J,R))
meanS=meanN+meanD
meanR=getPred(blockR[1,],sub,item,lag,I,J,R)

wN=rtnorm(RN,meanN[cond==0],rep(1,RN),s.crit[cbind(sub[cond==0]+1,resp[cond==0]+1,1)],s.crit[cbind(sub[cond==0]+1,resp[cond==0]+2,1)])

wS[notk[cond==1]]=rtnorm(Rnotk,meanS[notk],rep(1,RN),s.crit[cbind(sub[notk]+1,resp[notk]+1,1)],s.crit[cbind(sub[notk]+1,resp[notk]+2,1)])
wR[notk[cond==1]]=rtnorm(Rnotk,meanR[notk],rep(1,Rnotk),rep(-Inf,Rnotk),rep(0,Rnotk))

ind=isk[cond==1] & wR>0
wS[ind]=rnorm(sum(ind),meanS[cond==1][ind],1)
ind=isk[cond==1] & wR<0
wS[ind]=rtnorm(sum(ind),meanS[cond==1][ind],rep(1,sum(ind)),s.crit[cbind(sub[cond==1][ind]+1,K,1)],rep(Inf,sum(ind)))

ind=isk[cond==1] & wS>s.crit[cbind(sub[cond==1]+1,K,1)]
wR[ind]=rnorm(sum(ind),meanR[cond==1][ind],1)
ind=isk[cond==1] & wS<s.crit[cbind(sub[cond==1]+1,K,1)]
wR[ind]=rtnorm(sum(ind),meanR[cond==1][ind],rep(1,sum(ind)),rep(0,sum(ind)),rep(Inf,sum(ind)))


#Sample Blocks
allwN[cond==0]=wN
allwN[cond==1]=wS-meanD[cond==1]
tmp=sampleNorm(blockN[1,],allwN,sub,item,lag,I,J,R,nsubT,nitemT,10,.01,.01,met[1,1],met[1,2],1,0)
blockN[2,]=tmp[[1]]
b0[1,]=b0[1,]+tmp[[2]]

meanN=getPred(blockN[2,],sub,item,lag,I,J,R)
wD=wS-meanN[cond==1]

#sample muD
prop=blockD[1,mu]+rnorm(1,0,metD[mu])
oldLike= cond.muD(wD,blockD[1,mu],blockD[1,alpha],blockD[1,beta],blockD[1,theta],sub[cond==1],item[cond==1],lag[cond==1],sig2mu)
propLike=cond.muD(wD,prop          ,blockD[1,alpha],blockD[1,beta],blockD[1,theta],sub[cond==1],item[cond==1],lag[cond==1],sig2mu)
accept= rbinom(1,1,pmin(1,exp(propLike-oldLike)))
blockD[2,mu]=ifelse(accept,prop,blockD[1,mu])
b0D[mu]=b0D[mu]+accept

#sample alphaD
prop=blockD[1,alpha]+rnorm(I,0,metD[alpha])
oldLike= cond.alphaD(wD,blockD[2,mu],blockD[1,alpha],blockD[1,beta],blockD[1,theta],blockR[1,alpha],blockR[1,s2alpha],phi[1],sub[cond==1],item[cond==1],lag[cond==1],blockD[1,s2alpha])
propLike=cond.alphaD(wD,blockD[2,mu],prop             ,blockD[1,beta],blockD[1,theta],blockR[1,alpha],blockR[1,s2alpha],phi[1],sub[cond==1],item[cond==1],lag[cond==1],blockD[1,s2alpha])
accept= rbinom(I,1,pmin(1,exp(propLike-oldLike)))
blockD[2,alpha]=ifelse(accept,prop,blockD[1,alpha])
b0D[alpha]=b0D[alpha]+accept

#sample betaD
prop=blockD[1,beta]+rnorm(J,0,metD[beta])
oldLike= cond.betaD(wD,blockD[2,mu],blockD[2,alpha],blockD[1,beta],blockD[1,theta],blockR[1,beta],blockR[1,s2beta],phi[1],sub[cond==1],item[cond==1],lag[cond==1],blockD[1,s2beta])
propLike=cond.betaD(wD,blockD[2,mu],blockD[2,alpha],prop            ,blockD[1,theta],blockR[1,beta],blockR[1,s2beta],phi[1],sub[cond==1],item[cond==1],lag[cond==1],blockD[1,s2beta])
accept= rbinom(J,1,pmin(1,exp(propLike-oldLike)))
blockD[2,beta]=ifelse(accept,prop,blockD[1,beta])
b0D[beta]=b0D[beta]+accept

#sample thetaD
prop=blockD[1,theta]+rnorm(1,0,metD[theta])
oldLike=cond.thetaD (wD,blockD[2,mu],blockD[2,alpha],blockD[2,beta],blockD[1,theta],wR,blockR[1,mu],blockR[1,alpha],blockR[1,beta],sub[cond==1],item[cond==1],lag[cond==1],sig2mu)
propLike=cond.thetaD(wD,blockD[2,mu],blockD[2,alpha],blockD[2,beta],prop             ,wR,blockR[1,mu],blockR[1,alpha],blockR[1,beta],sub[cond==1],item[cond==1],lag[cond==1],sig2mu)
accept= rbinom(1,1,pmin(1,exp(propLike-oldLike)))
blockD[2,theta]=ifelse(accept,prop,blockD[1,theta])
b0D[theta]=b0D[theta]+accept


if(1==0){
#Decor MuD and AlphaD
shift=rnorm(1,0,metD[s2alpha])
muProp=blockD[m,mu]+shift
alphaProp=blockD[m,alpha]-shift
oldLike=-.5*sum(blockD[m,alpha]^2/blockD[m-1,s2alpha] + ((phi[m-1,1]*blockD[m,alpha])^2 - 2*phi[m-1,1]*blockD[m,alpha]*blockR[m-1,alpha])/blockR[m-1,s2alpha] + blockD[m,mu]^2/sig2mu)
propLike=-.5*sum(alphaProp^2/blockD[m-1,s2alpha] + ((phi[m-1,1]*alphaProp)^2 - 2*phi[m-1,1]*alphaProp*blockR[m-1,alpha])/blockR[m-1,s2alpha]+ muProp^2/sig2mu)
accept= rbinom(1,1,pmin(1,exp(propLike-oldLike)))
if(accept)
  {
    blockD[m,mu]=muProp
    blockD[m,alpha]=alphaProp
    b0D[s2alpha]=b0D[s2alpha]+1
  }

#Decor AlphaD and BetaD
shift=rnorm(1,0,metD[s2beta])
alphaProp=blockD[m,alpha]-shift
betaProp=blockD[m,beta]+shift
oldLike=-.5*sum(blockD[m,alpha]^2/blockD[m-1,s2alpha] + ((phi[m-1,1]*blockD[m,alpha])^2 - 2*phi[m-1,1]*blockD[m,alpha]*blockR[m-1,alpha])/blockR[m-1,s2alpha]) + -.5*sum(blockD[m,beta]^2/blockD[m-1,s2beta] + ((phi[m-1,2]*blockD[m,beta])^2 - 2*phi[m-1,2]*blockD[m,beta]*blockR[m-1,beta])/blockR[m-1,s2beta])
propLike=-.5*sum(alphaProp^2/blockD[m-1,s2alpha] + ((phi[m-1,1]*alphaProp)^2 - 2*phi[m-1,1]*alphaProp*blockR[m-1,alpha])/blockR[m-1,s2alpha]) + -.5*sum(betaProp^2/blockD[m-1,s2beta] + ((phi[m-1,2]*betaProp)^2 - 2*phi[m-1,2]*betaProp*blockR[m-1,beta])/blockR[m-1,s2beta])
accept= rbinom(1,1,pmin(1,exp(propLike-oldLike)))
if(accept)
  {
    blockD[m,beta]=betaProp
    blockD[m,alpha]=alphaProp
    b0D[s2beta]=b0D[s2beta]+1
  }
}

#SAMPLE VARIANCES
postB=sum(blockD[2,alpha]^2)/2 + b
blockD[2,s2alpha]=1/rgamma(1,shape=(a+I/2),scale=1/postB)
postB=sum(blockD[2,beta]^2)/2 + b
blockD[2,s2beta]=1/rgamma(1,shape=(a+J/2),scale=1/postB)

#Sample BlockR
tmp=sampleNormR(blockR[1,],phi[1,],blockD[2,],wR,sub[cond==1],item[cond==1],lag[cond==1],I,J,RS,nsubS,nitemS,sig2mu,.01,.01,met[2,1],met[2,2],1,1)
blockR[2,]=tmp[[1]]
b0[2,]=b0[2,]+tmp[[3]]

#Sample Phi

#########FIX 
tmpA=1/((sum(blockD[m,alpha]^2))/blockR[m,s2alpha] + 1/sig2mu)
tmpM=sum(blockR[m,alpha]*blockD[m,alpha])/blockR[m,s2alpha]
phi[m,1]=rnorm(1,tmpA*tmpM,sqrt(tmpA))

tmpA=1/((sum(blockD[m,beta]^2))/blockR[m,s2beta] + 1/sig2mu)
tmpM=sum(blockR[m,beta]*blockD[m,beta])/blockR[m,s2beta]
phi[m,2]=rnorm(1,tmpA*tmpM,sqrt(tmpA))


#####CAUTION!!!
#phi[2]=0

#Sample Criteria
PropCrit[,c(2,3,5,6)]=s.crit[,c(2,3,5,6),1]+rnorm(I*4,0,rep(met.crit,4))
violate=pmin(1,colSums((apply(PropCrit[,c(2,3,4,5,6)],1,diff))<0))
PropCrit[violate==1,]=s.crit[violate==1,,1]

likeOld=tapply(dpsdPosLogLike(R,I,J,K,resp,cond,sub,item,lag,blockN[2,],blockD[2,],blockR[2,],s.crit[,,1]),sub,sum)
likeProp=tapply(dpsdPosLogLike(R,I,J,K,resp,cond,sub,item,lag,blockN[2,],blockD[2,],blockR[2,],PropCrit),sub,sum)
accept=rbinom(I,1,pmin(1,(1-violate)*exp(likeProp-likeOld)))
b0.crit=b0.crit+accept
s.crit[accept==1,,2]=PropCrit[accept==1,]
s.crit[accept==0,,2]=s.crit[accept==0,,1]

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

#SET KEPT SAMPLES
blockN[1,]=blockN[2,]
blockD[1,]=blockD[2,]
blockR[1,]=blockR[2,]
s.crit[,,1]=s.crit[,,2]
phi[1]=phi[2]

if(m %in% keep)
  {
    ind=which(m==keep)
    save.blockN[ind,]=blockN[2,]
    save.blockD[ind,]=blockD[2,]
    save.blockR[ind,]=blockR[2,]
    save.crit[,,ind]=s.crit[,,2]
    save.phi[ind]=phi[2]
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
pb=txtProgressBar(min=1,max=NK,style=3,width=10)

D0=0
for(m in 1:NK)
  {
    D0=D0+sum(-2*dpsdPosLogLike(R,I,J,K,resp,cond,sub,item,lag,save.blockN[m,],save.blockD[m,],save.blockR[m,],save.crit[,,m]))
    setTxtProgressBar(pb, m)
  }
D0=D0/NK
Dhat=sum(-2*dpsdPosLogLike(R,I,J,K,resp,cond,sub,item,lag,estN,estD,estR,estCrit))
pD=D0-Dhat
DIC=pD+D0
}


u=new("dpsd")
u@mu=mu
u@alpha=alpha
u@beta=beta
u@s2alpha=s2alpha
u@s2beta=s2beta
u@theta=theta
u@estN=estN
u@estS=estD
u@estR=c(estR,colMeans(phi[keep,]))
u@estCrit=estCrit
u@blockN=save.blockN
u@blockS=save.blockD
u@blockR=cbind(save.blockR,save.phi)
u@s.crit=save.crit
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
