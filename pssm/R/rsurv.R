rsurv <-
function(cov2,dat,m,accumulate=TRUE,rescale=1) {
#dat is a data frame,cov1 and cov2 are numeric indicating the indecies of the covariates for progression and mortality
#m is the number of intervals, not that max(times)<m
#parms is a m*(m+1)/2 + m +length(cov1)+length(cov2) vector
dat[,c('tdeath')]<-rescale*dat[,c('tdeath')]
nscov2=m+1
ncov2=length(cov2)
nparms=m+ncov2;
ms<-dim(dat)[1]
na=(dat$cdeath==1)
f1=floor(dat$tdeath)+1
f2=dat$tdeath-(f1-1)
mask0=(matrix(1:m,ms,m,byrow=TRUE)==matrix(f1,ms,m))
mask1=(matrix(1:m,ms,m,byrow=TRUE)<matrix(f1,ms,m))
outfcn<-function(parms) {
if(ncov2>0) lam=matrix(parms[1:m],ms,m,byrow=TRUE)+matrix(as.matrix(dat[,cov2])%*%matrix(parms[nscov2:(m+ncov2)],ncov2,1),ms,m)
else lam=matrix(parms[1:m],ms,m,byrow=TRUE)
elam=exp(lam)
if (accumulate){
	ssv=rowSums(mask0*matrix(dat$cdeath,ms,m)*lam)-rowSums(mask1*elam+(mask0*elam)*matrix(f2,ms,m))
	ss=sum(ssv)
	#print(c(parms,ss))
return(ss)
	}
else {
elam1=exp(rowSums(mask0*matrix(dat$cdeath,ms,m)*lam))
ss=elam1*exp(-rowSums(mask1*elam+(mask0*elam)*matrix(f2,ms,m)))
return(log(ss))}
}
return(outfcn)
}
