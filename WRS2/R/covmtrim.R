covmtrim<-function(x,tr=.2,p=length(x),grp=c(1:p)){
#
#  Estimate the covariance matrix for the sample trimmed means corresponding
#  to the data in the R variable x,
#  which is assumed to be stored in list mode or a matrix.
# (x[[1]] contains the data for group 1, x[[2]] the data for group 2, etc.)
#  The function returns a p by p matrix of covariances, the diagonal
#  elements being equal to the squared standard error of the sample
#  trimmed means, where p is the number of groups to be included.
#  By default, all the groups in x are used, but a subset of
#  the groups can be used via grp.  For example, if
#  the goal is to estimate the covariances between the sample trimmed
#  means for groups 1, 2, and 5, use the command grp<-c(1,2,5)
#  before calling this function.
#
#  The default amount of trimming is 20%
#
#  Missing values (values stored as NA) are not allowed.
#
#  This function uses winvar from chapter 2.
#
if(is.list(x))x=matl(x)
x=elimna(x)
x=listm(x)
if(!is.list(x))stop("The data are not stored in list mode or a matrix.")
p<-length(grp)
pm1<-p-1
for (i in 1:pm1){
ip<-i+1
if(length(x[[grp[ip]]])!=length(x[[grp[i]]]))stop("The number of observations in each group must be equal")
}
n<-length(x[[grp[1]]])
h<-length(x[[grp[1]]])-2*floor(tr*length(x[[grp[1]]]))
covest<-matrix(0,p,p)
covest[1,1]<-(n-1)*winvar(x[[grp[1]]],tr)/(h*(h-1))
for (j in 2:p){
jk<-j-1
covest[j,j]<-(n-1)*winvar(x[[grp[j]]],tr)/(h*(h-1))
for (k in 1:jk){
covest[j,k]<-(n-1)*wincor(x[[grp[j]]],x[[grp[k]]],tr)$cov/(h*(h-1))
covest[k,j]<-covest[j,k]
}
}
covmtrim<-covest
covmtrim
}