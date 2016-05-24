MS.test<-function(Y,X,reps){
orderm<-matrix(nrow=nrow(Y),ncol=ncol(Y))
for(i in 1:ncol(Y)){
orderm[,i]<-rank(Y[,i]) 
}
S<-matrix(ncol=1,nrow=max(X))
for(i in 1 :max(X)){
S[i]<-sum(orderm[X==i])/reps
}
N<-dim(Y)[1]*dim(Y)[2]
MS<-(12/(max(X)*(N+ncol(Y)))*sum(S^2))-3*(N+ncol(Y))
p.val<-pchisq(MS,max(X)-1,lower.tail = FALSE)
res<-data.frame(df=max(X)-1,MS=MS,Chi.sq_p.val=p.val)
colnames(res)<-c("df","MS.test.stat","P(Chi.sq>MS)")
res
}
