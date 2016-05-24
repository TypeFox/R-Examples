.packageName='hbmem'
	
#Nice plotting params.
defpar=function(r,c)
{
  list(mfrow=c(r,c),pty='s',pch=19,mar=c(4,4,3,.2),mgp=c(2.0,.6,0),cex=1.15,pty='s')
}

#Truncated Normal
rtnorm=function(N,mu,sigma,a,b)
  {
    y=1:N*0
    .C("rtruncnorm",as.double(y),as.integer(N),as.double(mu),as.double(sigma),as.double(a),as.double(b),NAOK=TRUE,PACKAGE=.packageName)[[1]]
  }

#Extract predicted from block
getPred=function(block,cond,sub,item,lag,N,I,J,R)
  {
    pred=1:R
tmp=.C("getPred",as.double(pred),as.double(block),as.integer(cond),as.integer(sub),as.integer(item),as.double(lag),as.integer(N),as.integer(I),as.integer(J),as.integer(R),PACKAGE=.packageName)
    return(tmp[[1]])
  }

#Simulate Hierarchical Normal Data
normalSim=function(N=1,I=30,J=300,mu=0,s2a=.2,s2b=.2,muS2=0,s2aS2=0,s2bS2=0)
{
  R=I*J
  alpha=rnorm(I,0,sqrt(s2a))
  beta=rnorm(J,0,sqrt(s2b))
  alphaS2=rnorm(I,0,sqrt(s2aS2))
  betaS2=rnorm(J,0,sqrt(s2bS2))

cond=sample(0:(N-1),R,replace=TRUE)
subj=rep(0:(I-1),each=J)
item=NULL
for(i in 1:I)
item=c(item,sample(0:(J-1),J,replace=FALSE))

  lag=rep(0,R)
  resp=1:R
    for(r in 1:R)
      {
        mean=mu[cond[r]+1]+alpha[subj[r]+1]+beta[item[r]+1]
        sd=sqrt(exp(muS2+alphaS2[subj[r]+1]+betaS2[item[r]+1]))
        resp[r]=rnorm(1,mean,sd)
      }
out=as.data.frame(cbind(cond,subj,item,lag,resp))      
return(list(out,alpha,beta))
} 

#BLOCK MEANS WITH ONE VARIANCE AND CONDITION EFFECTS
sampleNorm = function(sample,y,cond,subj,item,lag,N,I,J,R,ncond,nsub,nitem,s2mu,s2a,s2b,meta,metb,sigma2,sampLag=TRUE,Hier=TRUE)
{
b0=c(0,0)
tmp=.C("sampleNormal", as.double(sample),as.double(y),as.integer(cond),as.integer(subj),as.integer(item),as.double(lag),as.integer(N),as.integer(I),as.integer(J),as.integer(R),as.integer(ncond),as.integer(nsub),as.integer(nitem),as.double(s2mu),as.double(s2a),as.double(s2b),as.double(meta),as.double(metb),as.integer(b0),as.double(sigma2),as.integer(sampLag),as.integer(Hier),PACKAGE=.packageName)

samp=tmp[[1]]
b0=tmp[[19]]
return(list(samp,b0))
}


sampleNormR = function(sample,phi,blockD,y,subj,item,lag,I,J,R,nsub,nitem,s2mu,s2a,s2b,meta,metb,sigma2,sampLag)
{
b0=c(0,0)
tmp=.C("sampleNormalR", as.double(sample),as.double(phi),as.double(blockD),as.double(y),as.integer(subj),as.integer(item),as.double(lag),as.integer(I),as.integer(J),as.integer(R),as.integer(nsub),as.integer(nitem),as.double(s2mu),as.double(s2a),as.double(s2b),as.double(meta),as.double(metb),as.integer(b0),as.double(sigma2),as.integer(sampLag),PACKAGE=.packageName)

samp=tmp[[1]]
p=tmp[[2]]
b0=tmp[[18]]
return(list(samp,p,b0))
}

#BLOCK MEANS WITH BLOCK VARIANCES
sampleNormb = function(sample,y,cond,subj,item,lag,N,I,J,R,ncond,nsub,nitem,s2mu,s2a,s2b,meta,metb,blockSigma2,sampLag=1,Hier=1)
{
b0=c(0,0)
tmp=.C("sampleNormalb", as.double(sample),as.double(y),as.integer(cond),as.integer(subj),as.integer(item),as.double(lag),as.integer(N),as.integer(I),as.integer(J),as.integer(R),as.integer(ncond),as.integer(nsub),as.integer(nitem),as.double(s2mu),as.double(s2a),as.double(s2b),as.double(meta),as.double(metb),as.integer(b0),as.double(blockSigma2),as.integer(sampLag),as.integer(Hier),NAOK=TRUE,PACKAGE=.packageName)

samp=tmp[[1]]
b0=tmp[[19]]
return(list(samp,b0))
}

#BLOCK VARIANCES
sampleSig2b = function(sample,y,cond,sub,item,lag,N,I,J,R,ncond,nsub,nitem,s2mu,s2a,s2b,met,blockMean,sampLag=1,Hier=1)
{
b0=rep(0,N+I+J+3)
tmp=.C("sampleSigma2b", as.double(sample),as.double(y),as.integer(cond),as.integer(sub),as.integer(item),as.double(lag),as.integer(N),as.integer(I),as.integer(J),as.integer(R),as.integer(ncond),as.integer(nsub),as.integer(nitem),as.double(s2mu),as.double(s2a),as.double(s2b),as.double(met),as.integer(b0),as.double(blockMean),as.integer(sampLag),as.integer(Hier),PACKAGE=.packageName)

samp=tmp[[1]]
b0=tmp[[18]]
return(list(samp,b0))
}


#ONE VARIANCE
sampleSig2=function(sig2,block,y,cond,sub,item,lag,N,ncond,I, J,a,b) 
return(.C("sampleSigma2",as.double(sig2),as.double(block),as.double(y),as.integer(cond),as.integer(sub),as.integer(item),as.double(lag),as.integer(N),as.integer(ncond),as.integer(I),as.integer(J),as.double(a),as.double(b),PACKAGE=.packageName)[[1]])


###########################
#NORMAL WITH POSITIVE MEAN#
###########################
nplike=function(x,mu,alpha,beta,theta,cond,sub,item,lag,sigma2)
{
mean=exp(mu[cond+1]+alpha[sub+1]+beta[item+1]+theta*lag)
return((mean^2 - 2*x*mean)/(-2*sigma2))
}

ld.mu=function(x,mu,alpha,beta,theta,cond,sub,item,lag,sig2mu,sigma2)
return(tapply(nplike(x,mu,alpha,beta,theta,cond,sub,item,lag,sigma2),cond,sum) - .5*(mu^2/sig2mu))

ld.alpha=function(x,mu,alpha,beta,theta,cond,sub,item,lag,sig2alpha,sigma2)
return(tapply(nplike(x,mu,alpha,beta,theta,cond,sub,item,lag,sigma2),sub,sum) - .5*(alpha^2/sig2alpha))

ld.beta=function(x,mu,alpha,beta,theta,cond,sub,item,lag,sig2beta,sigma2)
return(tapply(nplike(x,mu,alpha,beta,theta,cond,sub,item,lag,sigma2),item,sum) - .5*(beta^2/sig2beta))

ld.theta=function(x,mu,alpha,beta,theta,cond,sub,item,lag,sig2theta,sigma2)
return(sum(nplike(x,mu,alpha,beta,theta,cond,sub,item,lag,sigma2)) - .5*(theta^2/sig2theta))


samplePosNorm=function(sample,y,cond,sub,item,lag,N,I,J,R,sig2mu,a,b,met,sigma2,sampLag)
  {
    b0=rep(0,length(met))
    s.mu=sample[1:N]
    s.alpha=sample[(N+1):(N+I)]
    s.beta=sample[(N+I+1):(N+I+J)]
    s.s2alpha=sample[N+I+J+1]
    s.s2beta=sample[N+I+J+2]
    s.theta=sample[N+I+J+3]

  
#SAMPLE VARIANCES
postB=sum(s.alpha^2)/2 + b
s.s2alpha=1/rgamma(1,shape=(a+I/2),scale=1/postB)
postB=sum(s.beta^2)/2 + b
s.s2beta=1/rgamma(1,shape=(a+J/2),scale=1/postB)

#SAMPLE MU
prop=s.mu+rnorm(N,0,met[1:N])
oldLike=ld.mu(y,s.mu, s.alpha,s.beta,s.theta,cond,sub,item,lag,sig2mu,sigma2)
propLike=ld.mu(y,prop,s.alpha,s.beta,s.theta,cond,sub,item,lag,sig2mu,sigma2)
accept= rbinom(N,1,pmin(1,exp(propLike-oldLike)))
s.mu=ifelse(accept,prop,s.mu)
b0[1:N]=b0[1:N]+accept
#SAMPLE ALPHA
prop=s.alpha+rnorm(I,0,met[(N+1):(N+I)])
oldLike= ld.alpha(y,s.mu,s.alpha,s.beta,s.theta,cond,sub,item,lag,s.s2alpha,sigma2)
propLike=ld.alpha(y,s.mu,prop,   s.beta,s.theta,cond,sub,item,lag,s.s2alpha,sigma2)
accept= rbinom(I,1,pmin(1,exp(propLike-oldLike)))
s.alpha=ifelse(accept,prop,s.alpha)
b0[(N+1):(N+I)]=b0[(N+1):(N+I)]+accept
#SAMPLE BETA
prop=s.beta+rnorm(J,0,met[(N+I+1):(N+I+J)])
oldLike=ld.beta (y,s.mu,s.alpha,s.beta,s.theta,cond,sub,item,lag,s.s2beta,sigma2)
propLike=ld.beta(y,s.mu,s.alpha,prop  ,s.theta,cond,sub,item,lag,s.s2beta,sigma2)
accept= rbinom(J,1,pmin(1,exp(propLike-oldLike)))
s.beta=ifelse(accept,prop,s.beta)
b0[(N+I+1):(N+I+J)]=b0[(N+I+1):(N+I+J)] + accept

#DECORR
#mu and alpha
shift=rnorm(1,0,met[N+I+J+1])
muProp=s.mu+shift
aProp=s.alpha-shift
p=(sum(s.alpha^2)-sum(aProp^2))/(2*s.s2alpha) + (sum(s.mu^2)-sum(muProp^2))/(2*sig2mu)
accept= rbinom(1,1,pmin(1,exp(p)))
if(accept)
  {
    s.mu=muProp
    s.alpha=aProp
    b0[N+I+J+1]=b0[N+I+J+1]+1
  }

#alpha and beta
shift=rnorm(1,0,met[N+I+J+2])
aProp=s.alpha+shift
bProp=s.beta-shift
p=(sum(s.alpha^2)-sum(aProp^2))/(2*s.s2alpha) + (sum(s.beta^2)-sum(bProp^2))/(2*s.s2beta) 
accept= rbinom(1,1,pmin(1,exp(p)))
if(accept)
  {
    s.alpha=aProp
    s.beta=bProp
    b0[N+I+J+2]=b0[N+I+J+2]+1
  }

#SAMPLE THETA
if(sampLag){
prop=s.theta+rnorm(1,0,met[N+I+J+3])
oldLike=ld.theta (y,s.mu,s.alpha,s.beta,s.theta,cond,sub,item,lag,sig2mu,sigma2)
propLike=ld.theta(y,s.mu,s.alpha,s.beta,prop   ,cond,sub,item,lag,sig2mu,sigma2)
accept= rbinom(1,1,pmin(1,exp(propLike-oldLike)))
s.theta=ifelse(accept,prop,s.theta)
b0[N+I+J+3]=b0[N+I+J+3] + accept
}
    
    sample[1:N]=s.mu
    sample[(N+1):(N+I)]=s.alpha
    sample[(N+I+1):(N+I+J)]=s.beta
    sample[N+I+J+1]=s.s2alpha
    sample[N+I+J+2]=s.s2beta
    sample[N+I+J+3]=s.theta

    return(list(sample,b0))
}
