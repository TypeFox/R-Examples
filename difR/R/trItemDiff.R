trItemDiff<-function(prop, anchor = 1:nrow(prop)){
prop2<-prop
prop2[prop2>0.999]<-0.999
prop2[prop2<0.001]<-0.001
deltas<-4*qnorm(1-prop2)+13
m<-colMeans(deltas[anchor,])
sig<-cov(deltas[anchor,])
b<-(sig[2,2]-sig[1,1]+sqrt((sig[2,2]-sig[1,1])^2+4*sig[1,2]^2))/(2*sig[1,2])
a<-m[2]-b*m[1]
dist<-(b*deltas[,1]+a-deltas[,2])/sqrt(b^2+1)
pars<-rbind(c(a,b))
colnames(pars)<-c("a","b")
return(list(prop=prop,delta=deltas,pars=pars,dist=dist))}


