### R code from vignette source 'choplumpValidation.Rnw'

###################################################
### code chunk number 1: PreliminaryCalculations
###################################################
#source("H:\\main\\methods\\choplump\\r\\chop.functions.R")
library(choplump)
SEED<-1:200


###################################################
### code chunk number 2: choplumpValidation.Rnw:71-72
###################################################
chooseMatrix(5,2)


###################################################
### code chunk number 3: choplumpValidation.Rnw:76-78
###################################################
chopGeneral
choplumpGeneral


###################################################
### code chunk number 4: choplumpValidation.Rnw:84-98
###################################################
make.data<-function(N,M,SEED){
    set.seed(SEED)
    Z<-rep(0,N)
    Z[sample(1:N,N/2,replace=FALSE)]<-1
    test<-data.frame(W=c(rep(0,N-M),abs(rnorm(M))),Z=Z)
    return(test)
}

test<-make.data(10,6,SEED[1])
test
choplumpGeneral(test$W,test$Z,testfunc=testfunc.wilcox.ties.general)
cout<-choplump(W~Z,data=test,use.ranks=TRUE,exact=TRUE)
cout
cout$p.values


###################################################
### code chunk number 5: choplumpValidation.Rnw:104-109
###################################################
testfunc.DiM.general<-function(d){
    -TDiM(d$W,d$Z)
}
choplumpGeneral(test$W,test$Z,testfunc=testfunc.DiM.general)
choplump(W~Z,data=test,use.ranks=FALSE,exact=TRUE)$p.values


###################################################
### code chunk number 6: choplumpValidation.Rnw:121-140
###################################################
library(coin)
test<-make.data(20,12,SEED[1])
test
wilcox.manyzeros.exact(W=test$W,Z=test$Z)
test2<-data.frame(W=test$W,Z=as.factor(test$Z))
wilcox_test(W~Z,data=test2,distribution="exact",alternative="less")
wilcox_test(W~Z,data=test2,distribution="exact",alternative="greater")
wilcox_test(W~Z,data=test2,distribution="exact",alternative="two.sided")

test<-make.data(1000,12,SEED[2])
t0<-proc.time()
wilcox.manyzeros.exact(W=test$W,Z=test$Z)
## time for our algorithm
proc.time()-t0
test2<-data.frame(W=test$W,Z=as.factor(test$Z))
t1<-proc.time()
wilcox_test(W~Z,data=test2,distribution="exact",alternative="two.sided")
## time for coin algorithm
proc.time()-t1


###################################################
### code chunk number 7: choplumpValidation.Rnw:154-171
###################################################
test<-function(N,M,SEED,Use.Ranks=TRUE){
    out.approx<-out.exact<-rep(NA,length(SEED))
    for (i in 1:length(SEED)){
        set.seed(SEED[i])
        test<-make.data(N,M,SEED[i])
        out.approx[i]<-choplump(W~Z,data=test,use.ranks=Use.Ranks,method="approx")$p.values[1]
        out.exact[i]<-choplump(W~Z,data=test,use.ranks=Use.Ranks,method="exact",printNumCalcs=FALSE)$p.values[1]
    }
    out<-data.frame(plower.approx=out.approx,plower.exact=out.exact)
    return(out)
}
tout.ranks<-test(100,10,1:200,Use.Ranks=TRUE)
#tout.ranks
tout.noranks<-test(100,10,1:200,Use.Ranks=FALSE)
#tout.noranks
qrank<-quantile(tout.ranks[,1]-tout.ranks[,2],probs=c(.025,.975))
qnor<-quantile(tout.noranks[,1]-tout.noranks[,2],probs=c(.025,.975))


###################################################
### code chunk number 8: figureApproxVsExact
###################################################
#plot(tout.ranks[,1],tout.ranks[,2],xlim=c(0,1),ylim=c(0,1),xlab="Approximate p-value",ylab="Exact p-value",pch=1,cex=1.2)
#points(tout.noranks[,1],tout.noranks[,2],pch=15,cex=1.2)
#lines(c(0,1),c(0,1),lty=2,col="gray")

plot(c(tout.ranks[,2],tout.noranks[,2]),c(tout.noranks[,1]-tout.noranks[,2],tout.ranks[,1]-tout.ranks[,2]),xlim=c(0,1),type="n",
    ylab="(Approximate p-value) - (Exact p-value)",xlab="Exact p-value")
points(tout.ranks[,2],tout.ranks[,1]-tout.ranks[,2],pch=1,cex=1.0)
points(tout.noranks[,2],tout.noranks[,1]-tout.noranks[,2],pch=2,cex=1.0)

lines(c(-1,2),rep(qrank[1],2),lty=1,col="gray")
lines(c(-1,2),rep(qrank[2],2),lty=1,col="gray")
lines(c(-1,2),rep(qnor[1],2),lty=2,col="gray")
lines(c(-1,2),rep(qnor[2],2),lty=2,col="gray")


#plot(c(.5*tout.ranks[,2]+.5*tout.ranks[,1],.5*tout.noranks[,2]+.5*tout.noranks[,1]),c(tout.noranks[,1]-tout.noranks[,2],tout.ranks[,1]-tout.ranks[,2]),xlim=c(0,1),type="n",ylab="Approximate p value - Exact p value",xlab="Average of Approimate and Exact p value")
#points(.5*tout.ranks[,2]+.5*tout.ranks[,1],tout.ranks[,1]-tout.ranks[,2],pch=1,cex=1.2)
#points(.5*tout.noranks[,2]+.5*tout.noranks[,1],tout.noranks[,1]-tout.noranks[,2],pch=2,cex=1.2)
#qrank<-quantile(tout.ranks[,1]-tout.ranks[,2],probs=c(.025,.975))
#qnor<-quantile(tout.noranks[,1]-tout.noranks[,2],probs=c(.025,.975))
#lines(c(-1,2),rep(qrank[1],2),lty=1,col="gray")
#lines(c(-1,2),rep(qrank[2],2),lty=1,col="gray")
#lines(c(-1,2),rep(qnor[1],2),lty=2,col="gray")
#lines(c(-1,2),rep(qnor[2],2),lty=2,col="gray")




