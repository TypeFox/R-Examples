#### Help.Functions
library(quadprog)
library(MASS)
hspace<-function(beta,c)
{
 p<-length(beta)
 sgn.b<-sign(beta)
 #a<-sgn.b==0
 #sgn.b[a]<-1
 ord.b<-rank(abs(beta),ties.method="first")
 w<-c*(1:p-1)+(1-c)
 h.beta<-w[ord.b]*sgn.b 
 return(h.beta)
}


H.creater<-function(H,beta,t,c,epsilon)
{
 nrow.H<-nrow(H)
 p<-length(beta)
 ### active sets
  for(i in 1:nrow.H)
  {
    if(abs((sum(H[i,-1]*beta[-1])-t))>epsilon)H[i,]<-rep(0,p)#hspace(beta,c)
  }
  h.new<-t(hspace(beta[-1],c))
  ### vertical
  vert.check<-rep(0,nrow.H)
  for(i in 1:nrow.H)
  {
  if(sum(abs(H[i,-1]))!=0){vert.check[i]<-sum(h.new*H[i,-1])/(sqrt(sum(h.new^2))*sqrt(sum(H[i,-1]^2)))}
  }
  a<-(1-vert.check)<epsilon
  H[a,]<-rep(0,p)
  b<-rowSums(abs(H))!=0
  if(all(b)==FALSE){H.top<-NULL}
  if(any(b)==TRUE){H.top<-H[b,-1]}#matrix(,ncol=p)}
    
    #if(sum(H[b,])==0)nrow.H.top<-0
    #H<-H.top
    #if(all(a)==FALSE)
    #{
      H<-rbind(H.top,h.new)
    #}
  return(H) 
}


l2norm<-function(vec)
{
l2length<-sqrt(sum(vec^2))
if(l2length>0){vec.norm<-vec/l2length}
if(l2length==0){vec.norm<-vec}
return(list("length"=l2length,"normed"=vec.norm))
}



fs.oscar<-function(X,y,beta.alt, H.new, bvec, meq=0,familiy,epsilon=1e-8)
{
#if(family$family=="binomial")y<-as.factor(y)
p<-ncol(X)
beta.alt<-beta.alt#glm(y~X-1,family)$coef
beta.alt.l2<-l2norm(beta.alt[-1])$length
delta<-Inf
#X<-cbind(1,X)
while(delta>epsilon)
{
eta<-X%*%beta.alt
dev.mat<-family$mu.eta(eta)
mu<-family$linkinv(eta)
s.mat<-family$variance(mu)
W<-diag(c(dev.mat^2/s.mat))
Dmat<-t(X)%*%W%*%X
dvec<-t(X)%*%W%*%(eta+diag(c(1/dev.mat))%*%(y-mu))
beta.neu<-solve.QP(Dmat, dvec, -t(H.new), -bvec, meq=0, factorized=FALSE)$solution
delta<-l2norm(beta.neu-beta.alt)$length/beta.alt.l2
 
#cat(delta,"\n")
beta.alt<-beta.neu
beta.alt.l2<-l2norm(beta.alt)$length
}
return(beta.neu)
}




#### Active Set Estimate
oscar.as<-function(X,y,ml,t,c,epsilon=1e-8,family)
{
beta.alt<-ml#glm(y~X-1,family)$coef
H.new.eq.H.old<-FALSE
H.old<-matrix(0,1,ncol(X))
bvec<-rep(t,1)

#BETA<-NULL
while(H.new.eq.H.old==FALSE)
{
H.new<-cbind(0,H.creater(H.old,beta.alt,t,c,epsilon))#,epsilon=epsilon)
bvec<-c(rep(t,nrow(H.new)))
if(family$family!="gaussian")
{
beta<-fs.oscar(X, y,beta.alt, H.new, bvec, meq=0,family,epsilon)
}
if(family$family=="gaussian")
{
Dmat<-t(X)%*%X
dvec<-t(X)%*%y
beta<-solve.QP(Dmat, dvec, -t(H.new), -bvec, meq=0, factorized=FALSE)$solution
}
#if(all(dim(H.old)==dim(H.new))){H.new.eq.H.old<-all(H.old==H.new)}
if(l2norm(beta-beta.alt)$length/l2norm(beta.alt)$length<epsilon){H.new.eq.H.old<-TRUE}
H.old<-H.new                
#BETA<-cbind(BETA,beta)
#cat(max(abs(beta%*%t(H.new))),"\n")
beta.alt<-beta
#matplot(t(BETA),type="l")
}
return(beta)
}


glm.oscar<-function(y.org,X.org,family,t.seq,c,epsilon)
{
y<-y.org
p<-ncol(X.org)
X<-scale(X.org,center=TRUE,scale=TRUE)
X<-cbind(1,X)

control=glm.control(epsilon = 1e-8, maxit = 1000, trace = FALSE)
ml<-glm(y~X-1,family,control=control)
kq.hspace<-hspace(ml$coef[-1],c)
tmax<-sum(kq.hspace*ml$coef[-1])
t.seq<-tmax*t.seq#(tmax*0.01,tmax*0.99,length=t.length)
t.length<-length(t.seq)
#if(length(t.length)==1)t.seq<-c(t.length)
BETA<-matrix(0,t.length,p+1)
for(i in 1:t.length)
{
BETA[i,]<-oscar.as(X,y,ml$coef,t=t.seq[i],c=c,epsilon,family)
cat(i,"\n")
}
BETA.org<-BETA[,-1]%*%diag(1/sd(X.org))
Inter.org<-BETA[,1]-BETA.org%*%colMeans(X.org)
BETA.org<-cbind(Inter.org,BETA.org)
return(list(
"Beta.std"=BETA,
"Beta"=BETA.org,
"family"=family,
"t.seq"=t.seq
)
)
}




predict.glm.oscar<-function(BETA,X,y)
{
X<-cbind(1,X)
n<-nrow(X)
t.length<-nrow(BETA$Beta)
FIT<-matrix(0,n,t.length)
ETA<-matrix(0,n,t.length)
DEV<-matrix(0,n,t.length)
SUMDEV<-rep(0,t.length)
wt<-rep(1,n)
for(i in 1:t.length)
{
ETA[,i]<-X%*%BETA$Beta[i,]
FIT[,i]<-BETA$family$linkinv(ETA[,i])
DEV[,i]<-BETA$family$dev.resids(y, FIT[,i], wt)
SUMDEV[i]<-sum(DEV[,i])
}
return(list(
"eta"=ETA,
"fit"=FIT,
"dev"=DEV,
"sumdev"=SUMDEV
)
)
}


cv.glm.oscar<-function(X,y,family,splitM,c.seq,t.seq,epsilon)
{
nsplit<-ncol(splitM)
dev.all<-Inf
#dev<-rep(0,length(t.seq))
ncseq<-length(c.seq)
for(j in 1:ncseq)
{
dev<-rep(0,length(t.seq))
for(i in 1:nsplit)
{
X.train<-X[splitM[,i]==1,]
y.train<-y[splitM[,i]==1]
train<-glm.oscar(y.train,X.train,family=family,t.seq=t.seq,c=c.seq[j],epsilon=1e-8)
X.test<-X[splitM[,i]==0,]
y.test<-y[splitM[,i]==0]
test<-predict.glm.oscar(train,X.test,y.test)
dev<-dev+test$sumdev
}
cat(min(dev),"\n")
if(min(dev)<dev.all)
{
t.opt<-which.min(dev)
c.opt<-c.seq[j]
dev.all<-min(dev)

}

}
train<-glm.oscar(y,X,family=family,t.seq=t.seq[t.opt],c=c.opt,epsilon=1e-8)
return(list(
"model"=train,
"dev"=dev,
"c"=c.opt,
"t"=t.opt))
}
