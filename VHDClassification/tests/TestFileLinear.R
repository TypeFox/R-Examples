# TODO: Add comment
# 
# Author: robin
###############################################################################

library(VHDClassification)
############ Tests 2 classes linear and quadratic. What do we loose by searching a quadratic rule
############ when the true is linear
p=100; n=50 ; mu=array(0,c(p,2)); mu[1:10,1]=2 ;C=array(c(1,20),p)
xl=NULL; yl=NULL;
for (k in 1:2){xl=rbind(xl,t(array(C^(1/2),c(p,n))*(matrix(rnorm(p*n),nrow=p,ncol=n))+array(mu[,k],c(p,n))));
    yl=c(yl,array(k,n))}

x=NULL; y=NULL;
for (k in 1:2){
    x=rbind(x,t(array(C^(1/2),c(p,n))*(matrix(rnorm(p*n),nrow=p,ncol=n))+array(mu[,k],c(p,n))));
    y=c(y,array(k,n))    
}

#Learning
LearnedBinaryRule=learnBinaryRule(xl,yl,type='quadratic')
plotClassifRule(LearnedBinaryRule)
show(LearnedBinaryRule)
predict(LearnedBinaryRule,x)
LearnedBinaryRule=learnBinaryRule(xl,yl,type='linear')
predict(LearnedBinaryRule,x)
plotClassifRule(LearnedBinaryRule)
show(LearnedBinaryRule)
LearnedBinaryRule=learnBinaryRule(xl,yl,type='quadratic',procedure='UnivThresh')
plotClassifRule(LearnedBinaryRule)
show(LearnedBinaryRule)
predict(LearnedBinaryRule,x)
LearnedBinaryRule=learnBinaryRule(xl,yl,type='linear',procedure='UnivThresh')
predict(LearnedBinaryRule,x)
plotClassifRule(LearnedBinaryRule)
show(LearnedBinaryRule)
LearnedBinaryRule=learnBinaryRule(xl,yl,type='linear',procedure='FDRstudent')
predict(LearnedBinaryRule,x)
plotClassifRule(LearnedBinaryRule)
show(LearnedBinaryRule)
LearnedBinaryRule=learnBinaryRule(xl,yl,type='linear',procedure='Fisher')
predict(LearnedBinaryRule,x)
plotClassifRule(LearnedBinaryRule)
show(LearnedBinaryRule)
LearnedBinaryRule=learnBinaryRule(xl,yl,type='linear',procedure='FANThresh')
predict(LearnedBinaryRule,x)
plotClassifRule(LearnedBinaryRule)
show(LearnedBinaryRule)


############ Tests 3 classes linear
p=100; n=20 ; mu=array(0,c(p,4)); mu[1:10,1]=2 ;mu[11:20,2]=2;C=array(c(1,20),p)
mu[21:30,3]=2
x=NULL; y=NULL;
for (k in 1:4){
    x=rbind(x,t(array(C^(1/2),c(p,n))*(matrix(rnorm(p*n),nrow=p,ncol=n))+array(mu[,k],c(p,n))));
    y=c(y,array(k,n))}
#Learning
LearnedLinearPartitionWithLLR=learnPartitionWithLLR(x,y,procedure='FDRThresh')
LearnedQuadraticPartitionWithLLR=learnPartitionWithLLR(x,y,type='quadratic',procedure='FDRThresh')


plotClassifRule(LearnedLinearPartitionWithLLR)
#Testing Set
x=NULL; y=NULL;
for (k in 1:3){
    x=rbind(x,t(array(C^(1/2),c(p,n))*(matrix(rnorm(p*n),nrow=p,ncol=n))+array(mu[,k],c(p,n))));
    y=c(y,array(k,n))    
}
#Testing
myTestingdata=list(x=x,y=y)
LDAScore=mean((y!=predict(LearnedLinearPartitionWithLLR,myTestingdata$x))) ;
print(LDAScore)
QDAScore=mean((y!=predict(LearnedQuadraticPartitionWithLLR,myTestingdata$x))) ;
print(QDAScore)
