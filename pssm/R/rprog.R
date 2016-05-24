rprog <-
function(cov1,dat,m,accumulate=TRUE,scale=1,gradient=FALSE ,rescale=1) {
#dat is a data frame,cov1 and cov2 are numeric indicating the indecies of the covariates for progression and mortality
#m is the number of intervals, not that max(times)<m
#parms is a m*(m+1)/2 + m +length(cov1)+length(cov2) vector
dat[,c('tprog0','tprog1')]<-rescale*dat[,c('tprog0','tprog1')]
nscov1=m+1
ncov1=length(cov1)
nparms=m+ncov1;
ms<-dim(dat)[1]
na=!is.na(dat$tprog1)
m1=sum(na)
mask0<-matrix(1:m,ms,m,byrow=TRUE)<=matrix(floor(dat$tprog0)+1,ms,m)
if (m1>0) mask1<-(matrix(1:m,m1,m,byrow=TRUE)<=matrix(floor(dat$tprog1[na])+1,m1,m)) else mask1<-NULL

outfcn<-function(parms) {
if(ncov1>0) lam=matrix(parms[1:m],ms,m,byrow=TRUE)+matrix(as.matrix(dat[,cov1])%*%matrix(parms[nscov1:(m+ncov1)],ncov1,1),ms,m)
else lam=matrix(parms[1:m],ms,m,byrow=TRUE)
elam=exp(lam)
tmp0=    elam*(pmin(matrix(1:m,ms,m,byrow=TRUE),matrix(dat$tprog0,ms,m))-matrix(0:(m-1),ms,m,byrow=TRUE))
if (m1>0) tmp1=elam[na]*(pmin(matrix(1:m,m1,m,byrow=TRUE),matrix(dat$tprog1[na],m1,m))-matrix(0:(m-1),m1,m,byrow=TRUE))
ss=matrix(0,ms,1)
if (m1>0) ss[na]=-exp(-rowSums(tmp1*mask1))
ss= ss+exp(-rowSums(tmp0*mask0))
pdd=scale*log(ss)
if (accumulate) return(sum(pdd)) else return(pdd)
}
return(outfcn)
}
