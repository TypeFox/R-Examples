sampleNorm2 = function(sample,y,cond,subj,item,lag,N,I,J,R,ncond,nsub,nitem,s2mu,s2a,s2b,meta,metb,sigma2,sampLag=1,Hier=1)
{
b0=c(0,0)
tmp=.C("sampleNormal2", as.double(sample),as.double(y),as.integer(subj),as.integer(N),as.integer(I),as.integer(J),as.integer(R),as.integer(nsub),PACKAGE=.packageName)
samp=tmp[[1]]
return(list(samp))
}



dpsdSample2=function(dat,M=5000,keep=(M/10):M,getDIC=TRUE,Hier=TRUE,jump=.01)
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
if(!Hier)
  {
    blockN[,c(s2alphaN,s2betaN)]=2
    blockS[,c(s2alphaS,s2betaS)]=2
    blockR[,c(s2alphaS,s2betaS)]=.5
    cat("Non-Hierarchical Model!\n")
    cat("Using Prior Variances of 2.0 & .5 on Rec\n\n")
  }
wS=wR=rnorm(RS)
s.crit=array(dim=c(I,7,M))
s.crit[,1,]=-Inf
s.crit[,4,]=0
s.crit[,7,]=Inf
s.crit[,,1]=matrix(rep(c(-Inf,-1,-.5,0,.5,1,Inf),each=I),ncol=7)
b0=matrix(0,3,2)
met=matrix(0,3,2)
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


#Sample Blocks
tmp=sampleNorm2(blockN[m-1,],wN,cond[Scond==0],sub[Scond==0],item[Scond==0],lag[Scond==0],NN,I,J,RN,ncondN,nsubN,nitemN,100,.01,.01,met[1,1],met[1,2],1,1,Hier)
blockN[m,]=tmp[[1]]
#b0[1,]=b0[1,]+tmp[[2]]

tmp=sampleNorm2(blockS[m-1,],wS,cond[Scond==1],sub[Scond==1],item[Scond==1],lag[Scond==1],NS,I,J,RS,ncondS,nsubS,nitemS,100,.01,.01,met[2,1],met[2,2],1,1,Hier)
blockS[m,]=tmp[[1]]
#b0[2,]=b0[2,]+tmp[[2]]

tmp=sampleNorm2(blockR[m-1,],wR,cond[Scond==1],sub[Scond==1],item[Scond==1],lag[Scond==1],NS,I,J,RS,ncondS,nsubS,nitemS,100,.01,.01,met[3,1],met[3,2],1,1,Hier)
blockR[m,]=tmp[[1]]
#b0[3,]=b0[3,]+tmp[[2]]



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
    met.crit=met.crit+((b0.crit/m<.3)*-jump + (b0.crit/m>.5)*jump)
    met.crit[met.crit<jump]=jump
  }

if(m==keep[1])
  {
    cat("\n\nContinuing with following acceptance probabilities:\n")
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
    D0=D0+sum(-2*dpsdLogLike(R,NN,NS,I,J,K,resp,cond,Scond,sub,item,lag,blockN[m,],blockS[m,],blockR[m,],s.crit[,,m]))
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
