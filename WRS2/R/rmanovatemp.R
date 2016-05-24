rmanovatemp<-function(x,tr=.2,grp=c(1:length(x))){
#
#  A heteroscedastic one-way repeated measures ANOVA for trimmed means.
#
#  The data are assumed to be stored in $x$ which can
#  be either an n by J matrix, or an R variable having list mode.
#  If the data are stored in list mode,
#  length(x) is assumed to correspond to the total number of groups.
#  By default, the null hypothesis is that all group have a common mean.
#  To compare a subset of the groups, use grp to indicate which
#  groups are to be compared. For example, if you type the
#  command grp<-c(1,3,4), and then execute this function, groups
#  1, 3, and 4 will be compared with the remaining groups ignored.
#
if(!is.list(x) && !is.matrix(x))stop("Data must be stored in a matrix or in list mode.")
if(is.list(x)){
J<-length(grp)  # The number of groups to be compared
m1<-matrix(x[[grp[1]]],length(x[[grp[1]]]),1)
for(i in 2:J){     # Put the data into an n by J matrix
m2<-matrix(x[[grp[i]]],length(x[[i]]),1)
m1<-cbind(m1,m2)
}
}
if(is.matrix(x)){
if(length(grp)<ncol(x))m1<-as.matrix(x[,grp])
if(length(grp)>=ncol(x))m1<-as.matrix(x)
J<-ncol(x)
}
#
#  Raw data are now in the matrix m1
#
m2<-matrix(0,nrow(m1),ncol(m1))
xvec<-1
g<-floor(tr*nrow(m1))  #2g is the number of observations trimmed.
for(j in 1:ncol(m1)){  # Putting Winsorized values in m2
m2[,j]<-winval(m1[,j],tr)
xvec[j]<-mean(m1[,j],tr)
}
xbar<-mean(xvec)
qc<-(nrow(m1)-2*g)*sum((xvec-xbar)^2)
m3<-matrix(0,nrow(m1),ncol(m1))
m3<-sweep(m2,1,apply(m2,1,mean))  # Sweep out rows
m3<-sweep(m3,2,apply(m2,2,mean))  # Sweep out columns
m3<-m3+mean(m2)  # Grand Winsorized mean swept in
qe<-sum(m3^2)
test<-(qc/(qe/(nrow(m1)-2*g-1)))
#
#  Next, estimate the adjusted degrees of freedom
#
v<-winall(m1,tr=tr)$cov
vbar<-mean(v)
vbard<-mean(diag(v))
vbarj<-1
for(j in 1:J){
vbarj[j]<-mean(v[j,])
}
A<-J*J*(vbard-vbar)^2/(J-1)
B<-sum(v*v)-2*J*sum(vbarj^2)+J*J*vbar^2
ehat<-A/B
etil<-(nrow(m2)*(J-1)*ehat-2)/((J-1)*(nrow(m2)-1-(J-1)*ehat))
etil<-min(1.,etil)
df1<-(J-1)*etil
df2<-(J-1)*etil*(nrow(m2)-2*g-1)
siglevel<-1-pf(test,df1,df2)
list(test=test,df=c(df1,df2),siglevel=siglevel,tmeans=xvec,ehat=ehat,etil=etil)
}

