### R code from vignette source 'pp1.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: library
###################################################
library(polyapost)


###################################################
### code chunk number 2: R
###################################################
foo <- getRversion()
baz <- library(help = "polyapost")
baz <- baz$info[[1]]
baz <- baz[grep("Version", baz)]
baz <- sub("^Version: *", "", baz)
qux <- library(help = "rcdd")
qux <- qux$info[[1]]
qux <- qux[grep("Version", qux)]
qux <- sub("^Version: *", "", qux)


###################################################
### code chunk number 3: seed
###################################################
set.seed(313)


###################################################
### code chunk number 4: expp1
###################################################
ysamp<-c(0,1)
K<-10-2
polyap(ysamp,K)


###################################################
### code chunk number 5: expp2
###################################################
y<-rgamma(500,5)
mean(y)
samp<-sample(1:500,25)
ysamp<-y[samp]
mean(ysamp)
K<-500-25
simmns<-rep(0,10)
for(i in 1:10){simmns[i]<-mean(polyap(ysamp,K))}
round(simmns,digits=2)


###################################################
### code chunk number 6: expp3
###################################################
x<-sort(rgamma(500,10))
y<-rnorm(500,20 + 2*x,3)
cor(x,y)
mean(y)
mnx<-mean(x)
mnx


###################################################
### code chunk number 7: exp3.2
###################################################
samp<-sort(c(sample(1:250,8),sample(251:500,17)))
ysamp<-y[samp]
xsamp<-x[samp]
mean(ysamp)
mean(xsamp)
A1<-rbind(rep(1,25),c(rep(1,8),rep(0,17)))
b1<-c(1,0.5)
A2<-rbind(xsamp,-diag(25))
b2<-c(10.0,rep(0,25))
A3<-matrix(xsamp,1,25)
b3<-9.7


###################################################
### code chunk number 8: exp3.4
###################################################
eps<-0.001
initsol<-feasible(A1,A2,A3,b1,b2,b3,eps)
initsol[c(3,9,25)]


###################################################
### code chunk number 9: exp3.6
###################################################
burnin<-1
reps<-200001
out<-constrppmn(A1,A2,A3,b1,b2,b3,initsol,reps,ysamp,burnin)
mean(out[[1]])


###################################################
### code chunk number 10: figstrat
###################################################
plot(out[[1]][seq(1,200001,by=200)])
#plot(out[[1]])


###################################################
### code chunk number 11: exp3.7
###################################################
out[[2]] 
out[[3]] 


###################################################
### code chunk number 12: exp3.8
###################################################
A1<-rbind(rep(1,25))
A2<--diag(25)
b1<-1
b2<-rep(0,25)
initsol<-rep(0.04,25)
reps<-1000001
out<-constrppmn(A1,A2,NULL,b1,b2,NULL,initsol,reps,ysamp,burnin)
mean(out[[1]])
subseq<-seq(1,1000001,by=1000)
mean(out[[1]][subseq])
sqrt(var(out[[1]][subseq]))


###################################################
### code chunk number 13: figcprs
###################################################
plot(out[[1]][subseq])
#plot(out[[1]])


###################################################
### code chunk number 14: exp3.9
###################################################
K<-500-25
simmns<-rep(0,1000)
for(i in 1:1000){
       simmns[i]<-mean(polyap(ysamp,K))
}
mean(simmns)
sqrt(var(simmns))


###################################################
### code chunk number 15: figprs
###################################################
plot(simmns)


###################################################
### code chunk number 16: exp4
###################################################
A1<-rbind(rep(1,6),1:6)
A2<-rbind(c(2,5,7,1,10,8),diag(-1,6))
A3<-matrix(c(1,1,1,0,0,0),1,6)
b1<-c(1,3.5)
b2<-c(6,rep(0,6))
b3<-0.45
initsol<-rep(1/6,6)
out<-constrppprob(A1,A2,A3,b1,b2,b3,initsol,2000,5)
round(out,digits=5)


###################################################
### code chunk number 17: exp5
###################################################
ysamp<-c(1,2,3)
wts<-c(1,2,3)
wtpolyap(ysamp,wts,25)


###################################################
### code chunk number 18: pp1.Rnw:527-528
###################################################
smp<-c(20,10,10,25,15,20)


###################################################
### code chunk number 19: exphr
###################################################
 mxcst<-rbind(c(1,0,0,1,0,0),c(0,1,0,0,1,0),c(1,1,1,0,0,0))
 mncst<-c(5000/10000,3000/10000,6000/10000)
 out<-hitrun(smp, a2=mxcst, b2=mncst, nbatch=20, blen=1000)


###################################################
### code chunk number 20: exphr-means
###################################################
 round(colMeans(out$batch), digits=3)


###################################################
### code chunk number 21: exphr-mcse
###################################################
 round(apply(out$batch, 2, sd) / sqrt(out$nbatch), digits=3)


