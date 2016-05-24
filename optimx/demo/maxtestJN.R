# maxtestJN.R built on maxtestj.R
# check that optimx runs maximize for all methods

rm(list=ls())
library(optimx)

maxfn<-function(x) {
      	n<-length(x)
	ss<-seq(1,n)
	f<-10-(crossprod(x-ss))^2
	f<-as.numeric(f)
	return(f)
}


negmaxfn<-function(x) {
	f<-(-1)*maxfn(x)
	return(f)
}

x0<-rep(pi,4)
ans.mx<-optimx(x0,maxfn,control=list(maximize=TRUE,all.methods=TRUE,save.failures=TRUE,trace=TRUE))
print(ans.mx)
print(summary(ans.mx, order = value)[1,])

ans.mxn<-optimx(x0,negmaxfn,control=list(all.methods=TRUE,save.failures=TRUE,trace=TRUE))
print(ans.mxn)
print(summary(ans.mxn, order = value)[1,])


x00<-c(1,2,3,4)
# Test if things work when we provide the solution!
ans.mx0<-optimx(x0,maxfn,control=list(all.methods=TRUE,maximize=TRUE,save.failures=TRUE,trace=TRUE))
print(ans.mx0)
print(summary(ans.mx0, order = value)[1,])

ans.mx0n<-optimx(x0,negmaxfn,control=list(all.methods=TRUE,save.failures=TRUE,trace=TRUE))
print(ans.mx0n)
print(summary(ans.mx0n, order = value)[1,])



