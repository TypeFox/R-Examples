smean<-function(m,cop=3,MM=FALSE,op=1,outfun=outogk,cov.fun=rmba,MC=FALSE,STAND=FALSE,...){
#
# m is an n by p matrix
#
# Compute a multivariate skipped measure of location
#
# op=1:
# Eliminate outliers using a projection method
# If in addition, MC=T, a multi-core processor is used
# assuming your computer has multiple cores and the package
# multicore has been installed.
#
# cop=1 Donoho-Gasko median,
# cop=2 MCD,
# cop=3 marginal medians.
#  cop=4 uses MVE center
#  cop=5 uses TBS
#  cop=6 uses rmba (Olive's median ball algorithm)
#
# For each point
# consider the line between it and the center,
# project all points onto this line, and
# check for outliers using
#
# MM=F, a boxplot rule.
# MM=T, rule based on MAD and median
#
# Repeat this for all points. A point is declared
# an outlier if for any projection it is an outlier
# using a modification of the usual boxplot rule.
#
# op=2 use mgv (function outmgv) method to eliminate outliers
# an outlier if for any projection it is an outlier
# using a modification of the usual boxplot rule.
#
# op=3 use outlier method indicated by outfun
#
# Eliminate any outliers and compute means
#  using remaining data.
#
m<-elimna(m)
if(op==1){
if(!MC)temp<-outpro(m,plotit=FALSE,cop=cop,MM=MM,STAND=STAND)$keep
#if(MC)temp<-outproMC(m,plotit=FALSE,cop=cop,MM=MM,STAND=STAND)$keep
}
if(op==2)temp<-outmgv(m,plotit=FALSE,cov.fun=cov.fun)$keep
if(op==3)temp<-outfun(m,plotit=FALSE,...)$keep
val<-apply(m[temp,],2,mean)
val
}
