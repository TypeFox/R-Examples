johan<-function(cmat,vmean,vsqse,h,alpha=.05){
#
#  This function is used by other functions that come with this book,
#  and it can be used to test hypotheses not covered in the text.
#
#  The function performs Johansen's test of C mu = 0 for p independent groups,
#  where C is a k by p matrix of rank k and mu is a p by 1 matrix of
#  of unknown trimmed means.
#  The argument cmat contains the matrix C.
#  vmean is a vector of length p containing the p trimmed means
#  vsqe is a diagonal matrix containing the squared standard errors of the
#  the trimmed means in vmean.
#  h is a vector containing the effective sample sizes
#
yvec<-matrix(vmean,length(vmean),1)
if(!is.matrix(vsqse))vsqse<-diag(vsqse)
test<-cmat%*%vsqse%*%t(cmat)
invc<-solve(test)
test<-t(yvec)%*%t(cmat)%*%invc%*%cmat%*%yvec
R<-vsqse%*%t(cmat)%*%invc%*%cmat
A<-sum(diag((diag(R))^2/diag(h-1)))
df<-nrow(cmat)
crit<-qchisq(1-alpha,df)
crit<-crit+(crit/(2*df))*A*(1+3*crit/(df+2))
list(teststat=test[1],crit=crit[1])
}