johansp<-function(cmat,vmean,vsqse,h,J,K){
#
#  This function is used by other functions that come with this book,
#  and it can be used to test hypotheses not covered in the text.
#
#  The function performs Johansen's test of C mu = 0 for
#  a split-plot design where the first factor has independent groups,
#  while the second factor is within subjects,
#  C is a k by p matrix of rank k and mu is a p by 1 matrix of
#  of unknown trimmed means.
#  The argument cmat contains the matrix C.
#  vmean is a vector of length p containing the p trimmed means
#  vsqe is a block diagonal matrix, the jth block being the
#  estimated covariances among the trimmed means
#  in the jth level of factor A,
#  the trimmed means are in vmean,
#  h is a vector of length J containing the effective sample sizes for
#  the jth level of factor A.
#
p<-J*K
yvec<-matrix(vmean,length(vmean),1)
test<-cmat%*%vsqse%*%t(cmat)
invc<-solve(test)
test<-t(yvec)%*%t(cmat)%*%invc%*%cmat%*%yvec
temp<-0
klim<-1-K
kup<-0
for (j in 1:J){
klim<-klim+K
kup<-kup+K
Q<-matrix(0,p,p) #  create Q sub j
for (k in klim:kup)Q[k,k]<-1
mtem<-vsqse%*%t(cmat)%*%invc%*%cmat%*%Q
temp[j]<-(sum(diag(mtem%*%mtem))+(sum(diag(mtem)))^2)/(h[j]-1)
}
A<-.5*sum(temp)
df1<-nrow(cmat)
df2<-nrow(cmat)*(nrow(cmat)+2)/(3*A)
cval<-nrow(cmat)+2*A-6*A/(nrow(cmat)+2)
test<-test/cval
sig<-1-pf(test,df1,df2)
list(teststat=test[1],siglevel=sig)
}
