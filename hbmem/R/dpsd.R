setClass("dpsd",representation(
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
                               estR="numeric",
                               estCrit="matrix",
                               blockN="matrix",
                               blockS="matrix",
                               blockR="matrix",
                               s.crit="array",
                               pD="numeric",
                               DIC="numeric",
                               M="numeric",
                               keep="numeric",
                               b0="matrix",
                               b0Crit="numeric"))

setClass("dpsdSim",representation(
 			          Scond="numeric",                         
				  cond="numeric",
                                  subj="numeric",
                                  item="numeric",
                                  lag="numeric",
                                  resp="numeric",
                                  muN="numeric",
                                  muS="numeric",
                                  muR="numeric",
                                  alphaN="numeric",
                                  betaN="numeric",
                                  alphaS="numeric",
                                  betaS="numeric",
                                  alphaR="numeric",
                                  betaR="numeric"))

dpsdProbs=function(r,d,crit)
{
K=length(crit)+1
theta=diff(c(0,pnorm(crit,d,1),1))
ind=c(rep(0,K-1),1)
p=r*ind+(1-r)*theta
return(p)
}

dpsdSim=function(NN=2,NS=1,I=30,J=200,K=6,muN=c(-.7,-.5),s2aN=.2,s2bN=.2,muS=0,s2aS=.2,s2bS=.2,muR=qnorm(.25),s2aR=.2,s2bR=.2,crit=matrix(rep(c(-1.6,-.5,0,.5,1.6),each=I),ncol=(K-1)))
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

dpsdLogLike=function(R,NN,NS,I,JN,JS,K,dat,cond,Scond,sub,item,lag,blockN,blockS,blockR,crit)
  {
    l=1:R*0
.C("logLikeDpsd",as.double(l),as.integer(R),as.integer(NN),as.integer(NS),as.integer(I),as.integer(JN),as.integer(JS),as.integer(K),as.integer(dat),as.integer(cond),as.integer(Scond),as.integer(sub),as.integer(item),as.double(lag),as.double(blockN),as.double(blockS),as.double(blockR),as.double(as.vector(t(crit))),NAOK=TRUE,PACKAGE=.packageName)[[1]]
  }


dpsdSample=function(dat,M=5000,keep=(M/10):M,getDIC=TRUE,freeCrit=TRUE,Hier=TRUE,jump=.01)
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
NS=length(levels(as.factor(cond[Scond==1])))
I=length(levels(as.factor(sub)))
JN=length(levels(as.factor(item[Scond==0])))
JS=length(levels(as.factor(item[Scond==1])))
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
blockS=blockR=matrix(0,nrow=M,ncol=BS)
if(!Hier)
  {
    blockN[,c(s2alphaN,s2betaN)]=1
    blockS[,c(s2alphaS,s2betaS)]=1
    blockR[,c(s2alphaS,s2betaS)]=1
    cat("Non-Hierarchical Model!\n")
    cat("Using Prior Variances of 1\n\n")
  }
wS=wR=rnorm(RS)
s.crit=array(0,dim=c(I,K+1,M))
s.crit[,1,]=-Inf
s.crit[,K+1,]=Inf
fix.crit=floor(K/2)+1
var.crit=(2:K)[(2:K)!=fix.crit]
s.crit[,fix.crit,]=0
s.crit[,var.crit,1]=matrix(rep(seq(-.5,.5,length=length(var.crit)),each=I),ncol=length(var.crit))

b0=matrix(0,3,2)
met=matrix(.01,3,2)
met.crit=rep(.05,I)
b0.crit=rep(0,I)
PropCrit=s.crit[,,1]
if(I==1) PropCrit=t(as.matrix(s.crit[,,1]))
notk=Scond==1 & resp<(K-1)
Rnotk=sum(notk)
isk= Scond==1 & resp==(K-1) 
Risk=sum(isk)

pb=txtProgressBar(min=1,max=M,style=3,width=10)
cat("Starting MCMC\n")
#MCMC LOOP
for(m in 2:M){
#Sample latent data
wN=rtnorm(RN,getPred(blockN[m-1,],cond[Scond==0],sub[Scond==0],item[Scond==0],lag[Scond==0],NN,I,JN,RN),rep(1,RN),s.crit[cbind(sub[Scond==0]+1,resp[Scond==0]+1,m-1)],s.crit[cbind(sub[Scond==0]+1,resp[Scond==0]+2,m-1)])

wS[notk[Scond==1]]=rtnorm(Rnotk,getPred(blockS[m-1,],cond[notk],sub[notk],item[notk],lag[notk],NS,I,JS,Rnotk),rep(1,Rnotk),s.crit[cbind(sub[notk]+1,resp[notk]+1,m-1)],s.crit[cbind(sub[notk]+1,resp[notk]+2,m-1)])
wR[notk[Scond==1]]=rtnorm(Rnotk,getPred(blockR[m-1,],cond[notk],sub[notk],item[notk],lag[notk],NS,I,JS,Rnotk),rep(1,Rnotk),rep(-Inf,Rnotk),rep(0,Rnotk))

ind=isk[Scond==1] & wR>0
wS[ind]=rnorm(sum(ind),getPred(blockS[m-1,],cond[Scond==1][ind],sub[Scond==1][ind],item[Scond==1][ind],lag[Scond==1][ind],NS,I,JS,sum(ind)),1)
ind=isk[Scond==1] & wR<0
wS[ind]=rtnorm(sum(ind),getPred(blockS[m-1,],cond[Scond==1][ind],sub[Scond==1][ind],item[Scond==1][ind],lag[Scond==1][ind],NS,I,JS,sum(ind)),rep(1,sum(ind)),s.crit[cbind(sub[Scond==1][ind]+1,K,m-1)],rep(Inf,sum(ind)))

ind=isk[Scond==1] & wS>s.crit[cbind(sub[Scond==1]+1,K,m-1)]
wR[ind]=rnorm(sum(ind),getPred(blockR[m-1,],cond[Scond==1][ind],sub[Scond==1][ind],item[Scond==1][ind],lag[Scond==1][ind],NS,I,JS,sum(ind)),1)
ind=isk[Scond==1] & wS<s.crit[cbind(sub[Scond==1]+1,K,m-1)]
wR[ind]=rtnorm(sum(ind),getPred(blockR[m-1,],cond[Scond==1][ind],sub[Scond==1][ind],item[Scond==1][ind],lag[Scond==1][ind],NS,I,JS,sum(ind)),rep(1,sum(ind)),rep(0,sum(ind)),rep(Inf,sum(ind)))


#Sample Blocks
tmp=sampleNorm(blockN[m-1,],wN,cond[Scond==0],sub[Scond==0],item[Scond==0],lag[Scond==0],NN,I,JN,RN,ncondN,nsubN,nitemN,100,2,1,met[1,1],met[1,2],1,1,Hier)
blockN[m,]=tmp[[1]]
b0[1,]=b0[1,]+tmp[[2]]

tmp=sampleNorm(blockS[m-1,],wS,cond[Scond==1],sub[Scond==1],item[Scond==1],lag[Scond==1],NS,I,JS,RS,ncondS,nsubS,nitemS,100,2,1,met[2,1],met[2,2],1,1,Hier)
blockS[m,]=tmp[[1]]
b0[2,]=b0[2,]+tmp[[2]]

tmp=sampleNorm(blockR[m-1,],wR,cond[Scond==1],sub[Scond==1],item[Scond==1],lag[Scond==1],NS,I,JS,RS,ncondS,nsubS,nitemS,1,2,1,met[3,1],met[3,2],1,1,Hier)
blockR[m,]=tmp[[1]]
b0[3,]=b0[3,]+tmp[[2]]


  if(!freeCrit)
    {
      PropCrit[,var.crit]=s.crit[,var.crit,m-1]+rep(rnorm(length(var.crit),0,met.crit[1]),each=I)
      violate=pmin(1,sum(diff(PropCrit[1,2:K])<0))
      if(violate) PropCrit=t(as.matrix(s.crit[,,m-1]))

      likeOld=tapply(dpsdLogLike(R,NN,NS,I,JN,JS,K,resp,cond,Scond,sub,item,lag,blockN[m,],blockS[m,],blockR[m,],s.crit[,,m-1]),sub,sum)
      likeProp=tapply(dpsdLogLike(R,NN,NS,I,JN,JS,K,resp,cond,Scond,sub,item,lag,blockN[m,],blockS[m,],blockR[m,],PropCrit),sub,sum)
      accept=rbinom(I,1,pmin(1,(1-violate)*exp(likeProp-likeOld)))
      b0.crit=b0.crit+accept
      s.crit[,var.crit,m]=accept*PropCrit[,var.crit] + (1-accept)*s.crit[,var.crit,m-1]
    }

if(freeCrit)
  {
    PropCrit[,var.crit]=s.crit[,var.crit,m-1]+rnorm(length(var.crit)*I,0,rep(met.crit,length(var.crit)))
    violate=pmin(1,colSums((apply(PropCrit[,2:K],1,diff))<0))
    #PRIOR OF -7,7
    violate=pmin(1,violate + rowSums(PropCrit[,var.crit]<(-7)) + rowSums(PropCrit[,var.crit]>7))
 
   PropCrit[violate==1,]=s.crit[violate==1,,m-1]
    
    likeOld=tapply(dpsdLogLike(R,NN,NS,I,JN,JS,K,resp,cond,Scond,sub,item,lag,blockN[m,],blockS[m,],blockR[m,],s.crit[,,m-1]),sub,sum)
    likeProp=tapply(dpsdLogLike(R,NN,NS,I,JN,JS,K,resp,cond,Scond,sub,item,lag,blockN[m,],blockS[m,],blockR[m,],PropCrit),sub,sum)
    accept=rbinom(I,1,pmin(1,(1-violate)*exp(likeProp-likeOld)))

    b0.crit=b0.crit+accept
    s.crit[accept==1,,m]=PropCrit[accept==1,]
    s.crit[accept==0,,m]=s.crit[accept==0,,m-1]
  }
    
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
if(I==1) estCrit=as.matrix(rowMeans(s.crit[1,,keep]))

#GET DIC
DIC=pD=NA
if(getDIC)
{
print("Getting DIC")
pb=txtProgressBar(min=keep[1],max=tail(keep,1),style=3,width=10)

D0=0
for(m in keep)
  {
    D0=D0+sum(-2*dpsdLogLike(R,NN,NS,I,JN,JS,K,resp,cond,Scond,sub,item,lag,blockN[m,],blockS[m,],blockR[m,],s.crit[,,m]))
    setTxtProgressBar(pb, m)
  }
D0=D0/length(keep)
Dhat=sum(-2*dpsdLogLike(R,NN,NS,I,JN,JS,K,resp,cond,Scond,sub,item,lag,estN,estS,estR,estCrit))
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
