outbox<-function(x,mbox=FALSE,gval=NA,plotit=FALSE,STAND=FALSE){
#
# This function detects outliers using the
# boxplot rule, but unlike the R function boxplot,
# the ideal fourths are used to estimate the quartiles.
#
# Setting mbox=T results in using the modification
# of the boxplot rule suggested by Carling (2000).
#
x<-x[!is.na(x)] # Remove missing values
if(plotit)boxplot(x)
n<-length(x)
temp<-idealf(x)
if(mbox){
if(is.na(gval))gval<-(17.63*n-23.64)/(7.74*n-3.71)
cl<-median(x)-gval*(temp$qu-temp$ql)
cu<-median(x)+gval*(temp$qu-temp$ql)
}
if(!mbox){
if(is.na(gval))gval<-1.5
cl<-temp$ql-gval*(temp$qu-temp$ql)
cu<-temp$qu+gval*(temp$qu-temp$ql)
}
flag<-NA
outid<-NA
vec<-c(1:n)
for(i in 1:n){
flag[i]<-(x[i]< cl || x[i]> cu)
}
if(sum(flag)==0)outid<-NULL
if(sum(flag)>0)outid<-vec[flag]
keep<-vec[!flag]
outval<-x[flag]
n.out=sum(length(outid))
list(out.val=outval,out.id=outid,keep=keep,n=n,n.out=n.out,cl=cl,cu=cu)
}
