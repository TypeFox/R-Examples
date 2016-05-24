#cond.muD=function(x,mu,alpha,beta,theta,sub,item,lag,sig2mu)
#{
#mean=exp(mu+alpha[sub+1]+beta[item+1]+theta*lag)
#return(-.5*sum(mean^2 - 2*x*mean) + (mu^2/sig2mu))
#}

cond.rest.alphaD=function(wD,muD,alphaD,betaD,thetaD,wR,muR,alphaR,betaR,sub,item,lag,sig2alpha)
{
meanD=exp(muD+alphaD[sub+1]+betaD[item+1]+thetaD*lag)
likeD=-.5*(meanD^2 - 2*wD*meanD)
likeR=-.5*(alphaR[sub+1]^2 - 2*alphaR[sub+1]*(wR-muR-betaR[item+1]-thetaD*lag))
return(tapply(likeD,sub,sum) + tapply(likeR,sub,sum) - .5*(alphaD^2/sig2alpha)) 
}

cond.rest.betaD=function(wD,muD,alphaD,betaD,thetaD,wR,muR,alphaR,betaR,sub,item,lag,sig2beta)
{
meanD=exp(muD+alphaD[sub+1]+betaD[item+1]+thetaD*lag)
likeD=-.5*(meanD^2 - 2*wD*meanD)
likeR=-.5*(betaR[item+1]^2 - 2*betaR[item+1]*(wR-muR-alphaR[sub+1]-thetaD*lag))
return(tapply(likeD,item,sum) + tapply(likeR,item,sum) - .5*(betaD^2/sig2beta))
}

#cond.thetaD=function(wD,muD,alphaD,betaD,thetaD,wR,muR,alphaR,betaR,sub,item,lag,sig2theta)
#{
#meanD=exp(muD+alphaD[sub+1]+betaD[item+1]+thetaD*lag)
#likeD=-.5*(meanD^2 - 2*wD*meanD) 
#likeR=-.5*((thetaD*lag)^2 - 2*thetaD*lag*(wR-muR-alphaR[sub+1]-betaR[item+1]))
#return(sum(likeD)+sum(likeR)-.5*thetaD^2/sig2theta)
#}

dpsdRestSample=function(dat,M=5000,keep=(M/10):M,getDIC=TRUE,jump=.001)
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
blockN=blockD=blockR=matrix(0,nrow=M,ncol=B)
blockN[1,s2alpha]=blockN[1,s2beta]=blockD[1,s2alpha]=blockD[1,s2beta]=.01
phi=matrix(0,ncol=2,nrow=M)
wS=wR=rnorm(RS)
allwN=rnorm(R)
s.crit=array(dim=c(I,7,M))
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
meanN=getPred(blockN[m-1,],sub,item,lag,I,J,R)
meanD=exp(getPred(blockD[m-1,],sub,item,lag,I,J,R))
meanS=meanN+meanD
meanR=getPred(blockR[m-1,],sub,item,lag,I,J,R)

wN=rtnorm(RN,meanN[cond==0],rep(1,RN),s.crit[cbind(sub[cond==0]+1,resp[cond==0]+1,m-1)],s.crit[cbind(sub[cond==0]+1,resp[cond==0]+2,m-1)])

wS[notk[cond==1]]=rtnorm(Rnotk,meanS[notk],rep(1,RN),s.crit[cbind(sub[notk]+1,resp[notk]+1,m-1)],s.crit[cbind(sub[notk]+1,resp[notk]+2,m-1)])
wR[notk[cond==1]]=rtnorm(Rnotk,meanR[notk],rep(1,Rnotk),rep(-Inf,Rnotk),rep(0,Rnotk))

ind=isk[cond==1] & wR>0
wS[ind]=rnorm(sum(ind),meanS[cond==1][ind],1)
ind=isk[cond==1] & wR<0
wS[ind]=rtnorm(sum(ind),meanS[cond==1][ind],rep(1,sum(ind)),s.crit[cbind(sub[cond==1][ind]+1,K,m-1)],rep(Inf,sum(ind)))

ind=isk[cond==1] & wS>s.crit[cbind(sub[cond==1]+1,K,m-1)]
wR[ind]=rnorm(sum(ind),meanR[cond==1][ind],1)
ind=isk[cond==1] & wS<s.crit[cbind(sub[cond==1]+1,K,m-1)]
wR[ind]=rtnorm(sum(ind),meanR[cond==1][ind],rep(1,sum(ind)),rep(0,sum(ind)),rep(Inf,sum(ind)))

#Sample Blocks
allwN[cond==0]=wN
allwN[cond==1]=wS-meanD[cond==1]
tmp=sampleNorm(blockN[m-1,],allwN,sub,item,lag,I,J,R,nsubT,nitemT,10,.01,.01,met[1,1],met[1,2],1,0)
blockN[m,]=tmp[[1]]
b0[1,]=b0[1,]+tmp[[2]]

meanN=getPred(blockN[m,],sub,item,lag,I,J,R)
wD=wS-meanN[cond==1]
#sample muD
prop=blockD[m-1,mu]+rnorm(1,0,metD[mu])
oldLike= cond.muD(wD,blockD[m-1,mu],blockD[m-1,alpha],blockD[m-1,beta],blockD[m-1,theta],sub[cond==1],item[cond==1],lag[cond==1],sig2mu)
propLike=cond.muD(wD,prop          ,blockD[m-1,alpha],blockD[m-1,beta],blockD[m-1,theta],sub[cond==1],item[cond==1],lag[cond==1],sig2mu)
accept= rbinom(1,1,pmin(1,exp(propLike-oldLike)))
blockD[m,mu]=ifelse(accept,prop,blockD[m-1,mu])
b0D[mu]=b0D[mu]+accept

#sample alphaD
prop=blockD[m-1,alpha]+rnorm(I,0,metD[alpha])
oldLike= cond.rest.alphaD(wD,blockD[m,mu],blockD[m-1,alpha],blockD[m-1,beta],blockD[m-1,theta],wR,blockR[m-1,mu],blockR[m-1,alpha],blockR[m-1,beta],sub[cond==1],item[cond==1],lag[cond==1],blockD[m-1,s2alpha])
propLike=cond.rest.alphaD(wD,blockD[m,mu],prop             ,blockD[m-1,beta],blockD[m-1,theta],wR,blockR[m-1,mu],blockR[m-1,alpha],blockR[m-1,beta],sub[cond==1],item[cond==1],lag[cond==1],blockD[m-1,s2alpha])
accept= rbinom(I,1,pmin(1,exp(propLike-oldLike)))
blockD[m,alpha]=ifelse(accept,prop,blockD[m-1,alpha])
b0D[alpha]=b0D[alpha]+accept

#sample betaD
prop=blockD[m-1,beta]+rnorm(J,0,metD[beta])
oldLike= cond.rest.betaD(wD,blockD[m,mu],blockD[m,alpha],blockD[m-1,beta],blockD[m-1,theta],wR,blockR[m-1,mu],blockR[m-1,alpha],blockR[m-1,beta],sub[cond==1],item[cond==1],lag[cond==1],blockD[m-1,s2beta])
propLike=cond.rest.betaD(wD,blockD[m,mu],blockD[m,alpha],prop            ,blockD[m-1,theta],wR,blockR[m-1,mu],blockR[m-1,alpha],blockR[m-1,beta],sub[cond==1],item[cond==1],lag[cond==1],blockD[m-1,s2beta])
accept= rbinom(J,1,pmin(1,exp(propLike-oldLike)))
blockD[m,beta]=ifelse(accept,prop,blockD[m-1,beta])
b0D[beta]=b0D[beta]+accept

#sample thetaD
prop=blockD[m-1,theta]+rnorm(1,0,metD[theta])
oldLike=cond.thetaD (wD,blockD[m,mu],blockD[m,alpha],blockD[m,beta],blockD[m-1,theta],wR,blockR[m-1,mu],blockR[m-1,alpha],blockR[m-1,beta],sub[cond==1],item[cond==1],lag[cond==1],sig2mu)
propLike=cond.thetaD(wD,blockD[m,mu],blockD[m,alpha],blockD[m,beta],prop             ,wR,blockR[m-1,mu],blockR[m-1,alpha],blockR[m-1,beta],sub[cond==1],item[cond==1],lag[cond==1],sig2mu)
accept= rbinom(1,1,pmin(1,exp(propLike-oldLike)))
blockD[m,theta]=ifelse(accept,prop,blockD[m-1,theta])
b0D[theta]=b0D[theta]+accept


#SAMPLE VARIANCES
postB=sum(blockD[m,alpha]^2)/2 + b
blockD[m,s2alpha]=1/rgamma(1,shape=(a+I/2),scale=1/postB)
postB=sum(blockD[m,beta]^2)/2 + b
blockD[m,s2beta]=1/rgamma(1,shape=(a+J/2),scale=1/postB)


#Sample BlockR
blockR[m,theta]=blockD[m,theta]
blockR[m,c(s2alpha,s2beta)]=0
tmp=sum(wR-blockR[m-1,alpha[sub[cond==1]+1]]-blockR[m-1,beta[item[cond==1]+1]]-blockR[m,theta]*lag[cond==1])
tmpa=1/(RS + 1/sig2mu)
blockR[m,mu]=rnorm(1,tmp*tmpa,sqrt(tmpa))

#Sample Phis
tmpa=1/(sum(blockD[m,alpha[sub[cond==1]+1]]^2) + 1/sig2mu)
tmpM=sum(blockD[m,alpha[sub[cond==1]+1]]*(wR-blockR[m,mu]-blockR[m-1,beta[item[cond==1]+1]]-blockR[m,theta]*lag[cond==1]))
phi[m,1]=rnorm(1,tmpa*tmpM,sqrt(tmpa))
blockR[m,alpha]=phi[m,1]*blockD[m,alpha]

tmpa=1/(sum(blockD[m,beta[item[cond==1]+1]]^2) + 1/sig2mu)
tmpM=sum(blockD[m,beta[item[cond==1]+1]]*(wR-blockR[m,mu]-blockR[m,alpha[sub[cond==1]+1]]-blockR[m,theta]*lag[cond==1]))
phi[m,2]=rnorm(1,tmpa*tmpM,sqrt(tmpa))
blockR[m,beta] =phi[m,2]*blockD[m,beta]


#phi[m,]=c(2,4)

if(1==0){
#decor muD, alphaD, muR, alphaR
shift=rnorm(1,0,metD[s2alpha])
muDprop=blockD[m,mu]+shift
muRprop=blockR[m,mu]+shift
alphaProp=blockD[m,alpha]-shift
oldLike= -.5*sum(blockD[m,alpha]^2/blockD[m,s2alpha] + blockD[m,mu]^2/sig2mu + blockR[m,mu]^2/sig2mu)
propLike=-.5*sum(alphaProp^2/blockD[m,s2alpha] + muDprop^2/sig2mu + muRprop^2/sig2mu)
accept= rbinom(1,1,pmin(1,exp(propLike-oldLike)))
if(accept)
  {
    blockD[m,mu]=muDprop
    blockD[m,alpha]=alphaProp
    blockR[m,mu]=muRprop
    blockR[m,alpha]=phi[m,1]*alphaProp
    b0D[s2alpha]=b0D[s2alpha]+1
  }

#Decor AlphaD and BetaD
shift=rnorm(1,0,metD[s2beta])
alphaProp=blockD[m,alpha]-shift
betaProp=blockD[m,beta]+shift
p=(sum(blockD[m,alpha]^2)-sum(alphaProp^2))/(2*blockD[m-1,s2alpha]) + (sum(blockD[m,beta]^2)-sum(betaProp^2))/(2*blockD[m-1,s2beta])

accept= rbinom(1,1,pmin(1,exp(p)))
if(accept)
  {
    blockD[m,alpha]=alphaProp
    blockD[m,beta]=betaProp
    blockR[m,alpha]=phi[m,1]*alphaProp
    blockR[m,beta]=phi[m,2]*betaProp
    b0D[s2beta]=b0D[s2beta]+1
  }

}

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
u@mu=mu
u@alpha=alpha
u@beta=beta
u@s2alpha=s2alpha
u@s2beta=s2beta
u@theta=theta
u@estN=estN
u@estS=estD
u@estR=c(estR,colMeans(phi))
u@estCrit=estCrit
u@blockN=blockN[keep,]
u@blockS=blockD[keep,]
u@blockR=cbind(blockR[keep,],phi[keep,])
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
