gammaProbs=function(scale,shape,bounds)
{
cumulative=c(0,pgamma(bounds,shape=shape,scale=scale),1)
p.ij=diff(cumulative)
return(p.ij)
}

#Truncated Gamma
rtgamma=function(N,shape,scale,a,b)
  {
    y=1:N*0   
.C("rtruncgamma",as.double(y),as.integer(N),as.double(shape),as.double(scale),as.double(a),as.double(b),NAOK=TRUE,PACKAGE=.packageName)[[1]]
  }

gammaLogLike=function(R,NN,NS,I,J,K,dat,cond,Scond,subj,item,lag,blockN,blockS,crit,shape)
  {
    like=1:R*0
    
tmp=.C("logLikeGamma",as.double(like),as.integer(R),as.integer(NN),as.integer(NS),as.integer(I),as.integer(J),as.integer(K),as.integer(dat),as.integer(cond),as.integer(Scond),as.integer(subj),as.integer(item),as.double(lag),as.double(blockN),as.double(blockS),as.double(as.vector(t(crit))),as.double(shape),NAOK=TRUE,PACKAGE=.packageName)[[1]]
    return(tmp)
  }

gammaLikeLogLike=function(R,NN,NS,I,J,K,dat,cond,Scond,subj,item,lag,blockN,blockS,crit,shape)
  {
    like=1:R*0
    
tmp=.C("logLikeGammaLike",as.double(like),as.integer(R),as.integer(NN),as.integer(NS),as.integer(I),as.integer(J),as.integer(K),as.integer(dat),as.integer(cond),as.integer(Scond),as.integer(subj),as.integer(item),as.double(lag),as.double(blockN),as.double(blockS),as.double(as.vector(t(crit))),as.double(shape),NAOK=TRUE,PACKAGE=.packageName)[[1]]
    return(tmp)
  }

sampleGamma = function(sample,y,cond,subj,item,lag,N,I,J,R,ncond,nsub,nitem,s2mu,s2a,s2b,met,shape,sampLag,pos=FALSE)
{
  b0=rep(0,length(met))
if(pos==FALSE)  tmp=.C("sampleGamma", 
as.double(sample),as.double(y),as.integer(cond),as.integer(subj),as.integer(item),as.double(lag),as.integer(N),as.integer(I),as.integer(J),as.integer(R),as.integer(ncond),as.integer(nsub),as.integer(nitem),as.double(s2mu),as.double(s2a),as.double(s2b),as.double(met),as.integer(b0),as.double(shape),as.integer(sampLag),PACKAGE=.packageName)

if(pos==TRUE)  tmp=.C("samplePosGamma", 
as.double(sample),as.double(y),as.integer(cond),as.integer(subj),as.integer(item),as.double(lag),as.integer(N),as.integer(I),as.integer(J),as.integer(R),as.integer(ncond),as.integer(nsub),as.integer(nitem),as.double(s2mu),as.double(s2a),as.double(s2b),as.double(met),as.integer(b0),as.double(shape),as.integer(sampLag),PACKAGE=.packageName)
 
samp=tmp[[1]]
b0=tmp[[18]]
return(list(samp,b0))
}


gammaSim=function(NN=1,NS=2,I=30,J=200,K=6,muN=log(.65),s2aN=.2,s2bN=.2,muS=log(c(.8,1.2)),s2aS=.2,s2bS=.2,lagEffect=-.001,shape=2,crit=matrix(rep(c(.3,.6, 1, 1.2, 1.6),each=I),ncol=(K-1)))
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
        if(Scond[r]==0) p=gammaProbs(exp(muN[cond[r]+1]+alphaN[subj[r]+1]+betaN[item[r]+1]),2,crit[subj[r]+1,])
        if(Scond[r]==1) p=gammaProbs(exp(muS[cond[r]+1]+alphaS[subj[r]+1]+betaS[item[r]+1]+lag[r]*lagEffect),2,crit[subj[r]+1,])
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
    ret@muS2=0
    ret@alphaN=alphaN
    ret@alphaS=alphaS
    ret@alphaS2=0
    ret@betaN=betaN
    ret@betaS=betaS
    ret@betaS2=0
    ret@crit=crit
    return(ret)
  }



gammaLikeSim=function(NN=1,NS=1,I=30,J=200,K=6,muS=log(.5),s2aS=.2,s2bS=.2,lagEffect=0,shape=2,crit=matrix(rep(c(.75,.8,1,1.35,1.6),each=I),ncol=(K-1)))
  {
    muN=1
    #sanity checks:
    if(((J/2)/NN)%%1 !=0 | ((J/2)/NS)%%1 !=0) cat("Number of items must be divisible by number of conditions!\n")
 
    if(length(muN) !=NN) cat("Length of vector of new item means muN must match NN!\n")
    if(length(muS) !=NS) cat("Length of vector of studied item means muS must match NS!\n")

    R=I*J
    alphaS=rnorm(I,0,sqrt(s2aS))
    betaS=rnorm(J,0,sqrt(s2bS))
 
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
    scaleN=1

    dat=1:R
    for(r in 1:R)
      {
        scaleS=1+exp(muS[cond[r]+1]+alphaS[subj[r]+1]+betaS[item[r]+1]+lag[r]*lagEffect)
        if(sum(crit[subj[r]+1,]<(1/scaleS^2))) cat("Critera don't work!")
        
        if(Scond[r]==0)
          p=gammaProbs(scaleN,2,gammaLike2C(crit[subj[r]+1,],2,scaleN,scaleS))
            
        if(Scond[r]==1)
          p=gammaProbs(scaleS,2,gammaLike2C(crit[subj[r]+1,],2,scaleN,scaleS))

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
    ret@muS2=0
    ret@alphaN=0
    ret@alphaS=alphaS
    ret@alphaS2=0
    ret@betaN=0
    ret@betaS=betaS
    ret@betaS2=0
    ret@crit=crit
    return(ret)
  }


##############################################
##############################################
gammaSample=function(dat,M=10000,keep=(M/10):M,getDIC=TRUE,freeCrit=TRUE,shape=2,jump=.005)
{
Scond=dat$Scond
cond=dat$cond
sub=dat$sub
item=dat$item
lag=dat$lag-mean(dat$lag)
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
blockS=matrix(0,nrow=M,ncol=BS)
s.crit=array(dim=c(I,7,M))

blockN[1,c(s2alphaN,s2betaN)]=.1
blockS[1,c(s2alphaS,s2betaS)]=.1

s.crit[,1,]=0
s.crit[,4,]=1
s.crit[,7,]=Inf
s.crit[,,1]=matrix(rep(c(0,.3,.7,1,1.2,1.5,Inf),each=I),ncol=7)
b0S=rep(0,BS)
b0N=rep(0,BN)
metS=rep(.01,BS)
metN=rep(.01,BN)
met.crit=rep(.05,I)
b0.crit=rep(0,I)
PropCrit=s.crit[,,1]
if(I==1) PropCrit=t(as.matrix(s.crit[,,1]))

pb=txtProgressBar(min=1,max=M,style=3,width=10)
cat("Starting MCMC\n")
#MCMC LOOP
for(m in 2:M){
#Sample latent data
wN=rtgamma(RN,shape,exp(getPred(blockN[m-1,],cond[Scond==0],sub[Scond==0],item[Scond==0],lag[Scond==0],NN,I,J,RN)),s.crit[cbind(sub[Scond==0]+1,resp[Scond==0]+1,m-1)],s.crit[cbind(sub[Scond==0]+1,resp[Scond==0]+2,m-1)])

wS=rtgamma(RS,shape,exp(getPred(blockS[m-1,],cond[Scond==1],sub[Scond==1],item[Scond==1],lag[Scond==1],NS,I,J,RS)),s.crit[cbind(sub[Scond==1]+1,resp[Scond==1]+1,m-1)],s.crit[cbind(sub[Scond==1]+1,resp[Scond==1]+2,m-1)])
  
#Sample Blocks
tmp=sampleGamma(blockN[m-1,],wN,cond[Scond==0],sub[Scond==0],item[Scond==0],lag[Scond==0],NN,I,J,RN,ncondN,nsubN,nitemN,100,.01,.01,metN,shape,1)
blockN[m,]=tmp[[1]]
b0N=b0N+tmp[[2]]

tmp=sampleGamma(blockS[m-1,],wS,cond[Scond==1],sub[Scond==1],item[Scond==1],lag[Scond==1],NS,I,J,RS,ncondS,nsubS,nitemS,100,.01,.01,metS,shape,1)
blockS[m,]=tmp[[1]]
b0S=b0S+tmp[[2]]

#Sample Criteria
if(!freeCrit)
  {
    PropCrit[,c(2,3,5,6)]=s.crit[,c(2,3,5,6),m-1]+rep(rnorm(4,0,met.crit[1]),each=I)
    violate=pmin(1,sum(diff(PropCrit[1,1:6])<0))
    if(violate) PropCrit=s.crit[,,m-1]
    likeOld=sum(gammaLogLike(R,NN,NS,I,J,K,resp,cond,Scond,sub,item,lag,blockN[m,],blockS[m,],s.crit[,,m-1],shape))
    likeProp=sum(gammaLogLike(R,NN,NS,I,J,K,resp,cond,Scond,sub,item,lag,blockN[m,],blockS[m,],PropCrit,shape))
    accept=rbinom(1,1,min(1,(1-violate)*exp(likeProp-likeOld)))
    b0.crit=b0.crit+accept
    s.crit[,c(2,3,5,6),m]=accept*PropCrit[,c(2,3,5,6)] + (1-accept)*s.crit[,c(2,3,5,6),m-1]
  }

if(freeCrit)
{
  PropCrit[,c(2,3,5,6)]=s.crit[,c(2,3,5,6),m-1]+matrix(rnorm(I*4,0,rep(met.crit,4)),ncol=4)
  violate=pmin(1,colSums((apply(PropCrit[,1:6],1,diff))<0))
  PropCrit[violate==1,]=s.crit[violate==1,,m-1] 
  likeOld=tapply(gammaLogLike(R,NN,NS,I,J,K,resp,cond,Scond,sub,item,lag,blockN[m,],blockS[m,],s.crit[,,m-1],shape),sub,sum)
  likeProp=tapply(gammaLogLike(R,NN,NS,I,J,K,resp,cond,Scond,sub,item,lag,blockN[m,],blockS[m,],PropCrit,shape),sub,sum)
  accept=rbinom(I,1,pmin(1,(1-violate)*exp(likeProp-likeOld)))

  b0.crit=b0.crit+accept
  s.crit[accept==1,,m]=PropCrit[accept==1,]
  s.crit[accept==0,,m]=s.crit[accept==0,,m-1]
}

#AUTOTUNING
if(m>20 & m<keep[1] & m%%10==0)
  {
    metN=metN+         ((b0N/m<.3)*-jump + (b0N/m>.5)*jump)
    metS=metS+         ((b0S/m<.3)*-jump + (b0S/m>.5)*jump)
    met.crit=met.crit+((b0.crit/m<.3)*-jump + (b0.crit/m>.5)*jump)

    metN[metN<jump]=jump
    metS[metS<jump]=jump
    met.crit[met.crit<jump]=jump
  }

if(m==keep[1])
  {
    cat("\n\nFinal Tuning Parameters:\n")
    cat("New: ",round(metN,2),"\n",fill=TRUE)
    cat("Studied: ",round(metS,2),"\n",fill=TRUE)
    cat("Criteria: ",round(met.crit,2),"\n",fill=TRUE)


    cat("\n\nContinuing with following acceptance probabilities:\n")
    cat("New: ",round(b0N/m,2),"\n",fill=TRUE)
    cat("Studied: ",round(b0S/m,2),"\n",fill=TRUE)
    cat("Criteria: ",round(b0.crit/m,2),"\n",fill=TRUE)
    cat("If they are <.2 or >.6, adjust jump or keep\n")
  }
setTxtProgressBar(pb, m)
}

close(pb)

estN=colMeans(blockN[keep,])
estS=colMeans(blockS[keep,])
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
    D0=D0+sum(-2*gammaLogLike(R,NN,NS,I,J,K,resp,cond,Scond,sub,item,lag,blockN[m,],blockS[m,],s.crit[,,m],shape))
    setTxtProgressBar(pb, m)
  }
D0=D0/length(keep)
Dhat=sum(-2*gammaLogLike(R,NN,NS,I,J,K,resp,cond,Scond,sub,item,lag,estN,estS,estCrit,shape))
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
#u@estS2=0
u@estCrit=estCrit
u@blockN=blockN[keep,]
u@blockS=blockS[keep,]
#u@blockS2=0
u@s.crit=s.crit[,,keep]
u@pD=pD
u@DIC=DIC
u@M=M
u@keep=keep
u@b0=as.matrix(b0N/M)
u@b0S2=b0S/M
u@b0Crit=b0.crit/M

cat("Finished \n\n")
return(u)
}

#######################################################
#######################################################

gammaC2Like=function(crit,shape,thetaN,thetaS)
exp(crit*(1/thetaN - 1/thetaS) + shape*log(thetaN/thetaS)) 
  
gammaLike2C=function(like,shape,thetaN,thetaS)
(log(like) - shape*log(thetaN/thetaS))/(1/thetaN - 1/thetaS)

##############################################
##############################################
gammaLikeSample=function(dat,M=10000,keep=(M/10):M,getDIC=TRUE,shape=2,jump=.005)
{
Scond=dat$Scond
cond=dat$cond
sub=dat$sub
item=dat$item
lag=dat$lag-mean(dat$lag)
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

cat("\n CAUTION!!! ONLY WORKS FOR \n SINGLE CONDITION DESIGNS.\n\n")

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
blockS=matrix(0,nrow=M,ncol=BS)
blockS[1,muS]=rep(log(1.2),NS)
blockS[1,c(s2alphaS,s2betaS)]=0
s.crit=array(dim=c(I,7,M))
s.crit[,1,]=-Inf
s.crit[,7,]=Inf
s.crit[,,1]=matrix(rep(log(c(0,.65,.8,1,1.3,1.6,Inf)),each=I),ncol=7)
met.crit=rep(.05,I)
b0.crit=rep(0,I)
PropCrit=s.crit[,,1]

scale=like=likeProp=rep(0,R)
w=rep(0,R)
metS=rep(.01,BS)
b0S=rep(0,BS)
pb=txtProgressBar(min=1,max=M,style=3,width=10)
cat("Starting MCMC\n")


#MCMC LOOP
for(m in 2:M){
scaleS=1+exp(getPred(blockS[m-1,],cond,sub,item,lag,NS,I,J,R))
scale[Scond==0]=1-1/scaleS[Scond==0]
scale[Scond==1]=scaleS[Scond==1]-1
shift=2*log(scaleS)
lower=s.crit[cbind(sub+1,resp+1,m-1)] + shift
upper=s.crit[cbind(sub+1,resp+2,m-1)] + shift
w=rtgamma(R,shape,scale,lower,upper)-shift

#############
##SAMPLE MU##
#############

new=(w+2*log(scaleS))
likeOld=log(new) - new/scale -2*log(scale)
likeOld=tapply(likeOld,cond,sum) 

blockProp=blockS[m-1,]
blockProp[muS]=blockS[m-1,muS]+rnorm(1,0,metS[muS])
scaleS=1+exp(getPred(blockProp,cond,sub,item,lag,NS,I,J,R))
scale[Scond==0]=1-1/scaleS[Scond==0]
scale[Scond==1]=scaleS[Scond==1]-1
new=(w+2*log(scaleS))
likeProp=log(new) - new/scale -2*log(scale)
							
likeProp=tapply(likeProp,cond,sum) 
p=exp(likeProp-likeOld)
accept=rbinom(1,1,pmin(1,p))
if(sum(new<0)) {accept=0;cat(m,accept,blockProp[muS],"reject\n")}
b0S[muS]=b0S[muS]+accept
blockS[m,muS]=ifelse(accept,blockProp[muS],blockS[m-1,muS])

#blockS[m,muS]=blockS[1,muS]
#cat(blockS[m,muS],"\n")
blockS[m,2:BS]=blockS[m-1,2:BS]

###################
##Sample Criteria##
###################

if(1==1){
scaleS=1+exp(getPred(blockS[m,],cond,sub,item,lag,NS,I,J,R))
scale[Scond==0]=1-1/scaleS[Scond==0]
scale[Scond==1]=scaleS[Scond==1]-1

likeOld=log(pgamma(upper,shape,scale=scale)-pgamma(lower,shape,scale=scale))
likeOld=tapply(likeOld,sub,sum)

PropCrit[,2:6]=s.crit[,2:6,m-1]+matrix(rnorm(I*5,0,rep(met.crit,5)),ncol=5)
violate=pmin(1,colSums((apply(PropCrit[,2:6],1,diff))<0))
PropCrit[violate==1,]=s.crit[violate==1,,m-1] 
lower=PropCrit[cbind(sub+1,resp+1)] + shift
upper=PropCrit[cbind(sub+1,resp+2)] + shift
likeProp=log(pgamma(upper,shape,scale=scale)-pgamma(lower,shape,scale=scale))
likeProp=tapply(likeProp,sub,sum)

p=(1-violate)*exp(likeProp-likeOld)
accept=rbinom(I,1,pmin(1,p))

b0.crit=b0.crit+accept
s.crit[accept==1,,m]=PropCrit[accept==1,]
s.crit[accept==0,,m]=s.crit[accept==0,,m-1]
}

#min=as.vector(tapply(w,list(sub,resp),max)[,1:5])
#max=as.vector(tapply(w,list(sub,resp),min)[,2:6])
#s.crit[,2:6,m]=t(matrix(runif(I*(K-1),min,max),nrow=(K-1),byrow=T))

#s.crit[,,m]=s.crit[,,1]

#}

#par(mfrow=c(1,2))
#plot(exp(blockS[,1]),t='l')
#abline(h=.5)
#b0S[1]/M


#matplot(t(s.crit[1,,]),t='l')
#abline(h=t(s.crit[1,,1]))








##############
##AUTOTUNING##
if(m>20 & m<keep[1] & m%%10==0)
  {
#    metN=metN+         ((b0N/m<.3)*-jump + (b0N/m>.5)*jump)
    metS=metS+         ((b0S/m<.3)*-jump + (b0S/m>.5)*jump)
    met.crit=met.crit+((b0.crit/m<.3)*-jump + (b0.crit/m>.5)*jump)

#    metN[metN<jump]=jump
    metS[metS<jump]=jump
   met.crit[met.crit<jump]=jump
  }

if(m==keep[1])
  {
    cat("\n\nFinal Tuning Parameters:\n")
#    cat("New: ",round(metN,2),"\n",fill=TRUE)
    cat("Studied: ",round(metS,2),"\n",fill=TRUE)
    cat("Criteria: ",round(met.crit,2),"\n",fill=TRUE)


    cat("\n\nContinuing with following acceptance probabilities:\n")
#    cat("New: ",round(b0N/m,2),"\n",fill=TRUE)
    cat("Studied: ",round(b0S/m,2),"\n",fill=TRUE)
    cat("Criteria: ",round(b0.crit/m,2),"\n",fill=TRUE)
    cat("If they are <.2 or >.6, adjust jump or keep\n")
  }
setTxtProgressBar(pb, m)
}

close(pb)

estN=colMeans(blockN[keep,])
estS=colMeans(blockS[keep,])
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
    D0=D0+sum(-2*gammaLikeLogLike(R,NN,NS,I,J,K,resp,cond,Scond,sub,item,lag,blockN[m,],blockS[m,],s.crit[,,m],shape))
    setTxtProgressBar(pb, m)
  }
D0=D0/length(keep)
Dhat=sum(-2*gammaLikeLogLike(R,NN,NS,I,J,K,resp,cond,Scond,sub,item,lag,estN,estS,estCrit,shape))
pD=D0-Dhat
DIC=pD+D0
cat("\n\n")
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
#u@estN=estN
u@estS=estS
#u@estS2=0
u@estCrit=estCrit
#u@blockN=blockN[keep,]
u@blockS=blockS[keep,]
#u@blockS2=0
u@s.crit=s.crit[,,keep]
u@pD=pD
u@DIC=DIC
u@M=M
u@keep=keep
u@b0=as.matrix(b0S/M)
#u@b0S2=b0S/M
u@b0Crit=b0.crit/M

cat("Finished \n\n")
return(u)
}


