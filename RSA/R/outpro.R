#
# outpro<-function(m,gval=NA,center=NA,plotit=TRUE,op=TRUE,MM=FALSE,cop=3,
# xlab="VAR 1",ylab="VAR 2",STAND=TRUE,tr=.2,q=.5,pr=TRUE,...){
# #
# # Detect outliers using a modification of the
# # Stahel-Donoho  projection method.
# #
# # Determine center of data cloud, for each point,
# # connect it with center, project points onto this line
# # and use distances between projected points to detect
# # outliers. A boxplot method is used on the
# # projected distances.
# #
# # plotit=TRUE creates a scatterplot when working with
# # bivariate data.
# #
# # op=T
# # means the .5 depth contour is plotted
# # based on data with outliers removed.
# #
# # op=F
# # means .5 depth contour is plotted without removing outliers.
# #
# #  MM=F  Use interquatile range when checking for outliers
# #  MM=T  uses MAD.
# #
# #  If value for center is not specified,
# #  there are four options for computing the center of the
# #  cloud of points when computing projections:
# #
# #  cop=2 uses MCD center
# #  cop=3 uses median of the marginal distributions.
# #  cop=4 uses MVE center
# #  cop=5 uses TBS
# #  cop=6 uses rmba (Olive's median ball algorithm)#  cop=7 uses the spatial (L1) median
# #
# #  args q and tr having are not used by this function. They are included to deal
# #  with situations where smoothers have optional arguments for q and tr
# #
# #  When using cop=2, 3 or 4, default critical value for outliers
# #  is square root of the .975 quantile of a
# #  chi-squared distribution with p degrees
# #  of freedom.
# #
# #  STAND=T means that marginal distributions are standardized before
# #  checking for outliers.
# #
# #  Donoho-Gasko (Tukey) median is marked with a cross, +.
# #
# m<-as.matrix(m)
# if(pr){
# if(!STAND){
# if(ncol(m)>1)print("STAND=FALSE. If measures are on different scales, might want to use STAND=TRUE")
# }}
# library(MASS)
# m=elimna(m)
# m<-as.matrix(m)
# nv=nrow(m)
# if(ncol(m)==1){
# dis<-(m-median(m,na.rm=TRUE))^2/mad(m,na.rm=TRUE)^2
# dis<-sqrt(dis)
# dis[is.na(dis)]=0
# crit<-sqrt(qchisq(.975,1))
# chk<-ifelse(dis>crit,1,0)
# vec<-c(1:nrow(m))
# outid<-vec[chk==1]
# keep<-vec[chk==0]
# }
# if(ncol(m)>1){
# if(STAND)m=standm(m,est=median,scat=mad)
# if(is.na(gval) && cop==1)gval<-sqrt(qchisq(.95,ncol(m)))
# if(is.na(gval) && cop!=1)gval<-sqrt(qchisq(.975,ncol(m)))
# if(cop==1 && is.na(center[1])){
# if(ncol(m)>2)center<-dmean(m,tr=.5,cop=1)
# if(ncol(m)==2){
# tempd<-NA
# for(i in 1:nrow(m))
# tempd[i]<-depth(m[i,1],m[i,2],m)
# mdep<-max(tempd)
# flag<-(tempd==mdep)
# if(sum(flag)==1)center<-m[flag,]
# if(sum(flag)>1)center<-apply(m[flag,],2,mean)
# }}
# if(cop==2 && is.na(center[1])){
# center<-cov.mcd(m)$center
# }
# if(cop==4 && is.na(center[1])){
# center<-cov.mve(m)$center
# }
# if(cop==3 && is.na(center[1])){
# center<-apply(m,2,median)
# }
# if(cop==5 && is.na(center[1])){
# center<-tbs(m)$center
# }
# if(cop==6 && is.na(center[1])){
# center<-rmba(m)$center
# }
# if(cop==7 && is.na(center[1])){
# center<-spat(m)
# }
# flag<-rep(0, nrow(m))
# outid <- NA
# vec <- c(1:nrow(m))
# for (i in 1:nrow(m)){
# B<-m[i,]-center
# dis<-NA
# BB<-B^2
# bot<-sum(BB)
# if(bot!=0){
# for (j in 1:nrow(m)){
# A<-m[j,]-center
# temp<-sum(A*B)*B/bot
# dis[j]<-sqrt(sum(temp^2))
# }
# temp<-idealf(dis)
# if(!MM)cu<-median(dis)+gval*(temp$qu-temp$ql)
# if(MM)cu<-median(dis)+gval*mad(dis)
# outid<-NA
# temp2<-(dis> cu)
# flag[temp2]<-1
# }}
# if(sum(flag) == 0) outid <- NA
# if(sum(flag) > 0)flag<-(flag==1)
# outid <- vec[flag]
# idv<-c(1:nrow(m))
# keep<-idv[!flag]
# if(ncol(m)==2){
# if(plotit){
# plot(m[,1],m[,2],type="n",xlab=xlab,ylab=ylab)
# points(m[keep,1],m[keep,2],pch="*")
# if(length(outid)>0)points(m[outid,1],m[outid,2],pch="o")
# if(op){
# tempd<-NA
# keep<-keep[!is.na(keep)]
# mm<-m[keep,]
# for(i in 1:nrow(mm))tempd[i]<-depth(mm[i,1],mm[i,2],mm)
# mdep<-max(tempd)
# flag<-(tempd==mdep)
# if(sum(flag)==1)center<-mm[flag,]
# if(sum(flag)>1)center<-apply(mm[flag,],2,mean)
# m<-mm
# }
# points(center[1],center[2],pch="+")
# x<-m
# temp<-fdepth(m,plotit=FALSE)
# flag<-(temp>=median(temp))
# xx<-x[flag,]
# xord<-order(xx[,1])
# xx<-xx[xord,]
# temp<-chull(xx)
# #xord<-order(xx[,1])
# #xx<-xx[xord,]
# #temp<-chull(xx)
# lines(xx[temp,])
# lines(xx[c(temp[1],temp[length(temp)]),])
# }}}
# list(n=nv,n.out=length(outid),out.id=outid,keep=keep)
# }
#
#
# elimna<-function(m){
# #
# # remove any rows of data having missing values
# #
# if(is.null(dim(m)))m<-as.matrix(m)
# ikeep<-c(1:nrow(m))
# for(i in 1:nrow(m))if(sum(is.na(m[i,])>=1))ikeep[i]<-0
# elimna<-m[ikeep[ikeep>=1],]
# elimna
# }
#
# standm<-function(x,locfun=lloc,est=mean,scat=var,...){
# # standardize a matrix x
# #
# x=elimna(x)
# x=as.matrix(x)
# m1=lloc(x,est=est)
# v1=apply(x,2,scat)
# p=ncol(x)
# for(j in 1:p)x[,j]=(x[,j]-m1[j])/sqrt(v1[j])
# x
# }
#
# lloc<-function(x,est=tmean,...){
# if(is.data.frame(x)){
# x=as.matrix(x)
# x=apply(x,2,as.numeric) # earlier versions of R require this command
# }
# if(!is.list(x))val<-est(x,...)
# if(is.list(x))val=lapply(x,est)
# if(is.matrix(x))val<-apply(x,2,est,...)
# val
# }
#
#
# idealf<-function(x,na.rm=FALSE){
# #
# # Compute the ideal fourths for data in x
# #
# if(na.rm)x<-x[!is.na(x)]
# j<-floor(length(x)/4 + 5/12)
# y<-sort(x)
# g<-(length(x)/4)-j+(5/12)
# ql<-(1-g)*y[j]+g*y[j+1]
# k<-length(x)-j+1
# qu<-(1-g)*y[k]+g*y[k-1]
# list(ql=ql,qu=qu)
# }