out<-function(x,cov.fun=cov.mve,SEED=TRUE,xlab="X",ylab="Y",qval=.975,
crit=NULL,plotit=FALSE,...){
#
#  Search for outliers using robust measures of location and scatter,
#  which are used to compute robust analogs of Mahalanobis distance.
#
#  x is an n by p matrix or a vector of data.
#
#  The function returns the values flagged as an outlier plus
#  the (row) number where the data point is stored.
#  If x is a vector, out.id=4 indicates that the fourth observation
#  is an outlier and outval=123 indicates that 123 is the value.
#  If x is a matrix, out.id=4 indicates that the fourth row of
#  the matrix is an outlier and outval reports the corresponding
#  values.
#
#  The function also returns the distance of the
#  points identified as outliers
#  in the variable dis.
#
#  For bivariate data, if plotit=TRUE, plot points and circle outliers.
#
#  cov.fun determines how the measure of scatter is estimated.
#  Possible choices are
#  cov.mve (the MVE estimate)
#  cov.mcd (the MCD estimate)
#  covmba2 (the MBA or median ball algorithm)
#  rmba  (an adjustment of MBA suggested by D. Olive)
#  cov.roc (Rocke's TBS estimator)
#
#  plotit=FALSE used to avoid problems when other functions in WRS call 
#  this function
#

if(SEED)set.seed(12)
if(is.data.frame(x))x=as.matrix(x)
if(is.list(x))stop("Data cannot be stored in list mode")
nrem=nrow(as.matrix(x))
if(!is.matrix(x)){
dis<-(x-median(x,na.rm=TRUE))^2/mad(x,na.rm=TRUE)^2
if(is.null(crit))crit<-sqrt(qchisq(.975,1))
vec<-c(1:length(x))
}
if(is.matrix(x)){
mve<-cov.fun(elimna(x))
dis<-mahalanobis(x,mve$center,mve$cov)
if(is.null(crit))crit<-sqrt(qchisq(.975,ncol(x)))
vec<-c(1:nrow(x))
}
dis[is.na(dis)]=0
dis<-sqrt(dis)
chk<-ifelse(dis>crit,1,0)
id<-vec[chk==1]
keep<-vec[chk==0]
if(is.matrix(x)){
if(ncol(x)==2 && plotit){
plot(x[,1],x[,2],xlab=xlab,ylab=ylab,type="n")
flag<-rep(T,nrow(x))
flag[id]<-F
points(x[flag,1],x[flag,2])
if(sum(!flag)>0)points(x[!flag,1],x[!flag,2],pch="*")
}}
if(!is.matrix(x))outval<-x[id]
if(is.matrix(x))outval<-x[id,]
n=nrow(as.matrix(x))
n.out=length(id)
list(n=n,n.out=n.out,out.val=outval,out.id=id,keep=keep,dis=dis,crit=crit)
}
