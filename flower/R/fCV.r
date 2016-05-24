fCV<-function(dd,pop){
pop=as.factor(pop)
m=length(colnames(dd))
a=length(levels(pop))
popdd=array(0,dim=c(m,a))
colnames(popdd)=levels(pop)
rownames(popdd)=colnames(dd)
for (i in 1:m){
  popdd[i,]=tapply(dd[,i],pop,sum)
}
aij=array(0,dim=c(a,1))
rownames(aij)=levels(pop)
bij=array(0,dim=c(a,1))
rownames(bij)=levels(pop)
for (i in 1:a){
xtpdd=unique(popdd[,i])
aij[i,]=sd(xtpdd)
bij[i,]=mean(xtpdd)
}
cv=aij/bij
colnames(cv)=c("cv")
list(CV.x=popdd,CV.sd=aij,CV.mean=bij,CV=cv)
}