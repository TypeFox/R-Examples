# require library lme4

library(lme4)

# DIF: LIKELIHOOD-RATIO TEST

LRT<-function(data,member){
N<-nrow(data)
C<-ncol(data)
y<-NULL
for (i in 1:ncol(data)) y<-c(y,as.numeric(data[,i]))

person<-1:N
pp<-rep(person,C)
pp<-as.factor(pp)

member<-rep(member,C)
items<-rep(1,N)
for (i in 2:C) items<-c(items,rep(i,N))
items<-as.factor(items)


res<-NULL
nodif<-glmer(y ~items + (1|pp) -1 + member,family=binomial, REML=FALSE)

for (item in 1:C){
X<-rep(0,N*C)
mi<-(item-1)*N+1
ma<-item*N
X[mi:ma]<-1
dif<-glmer(y ~items + (1|pp) -1 + member + member:X,family=binomial, REML=FALSE)

res[item]<-deviance(nodif)-deviance(dif)
}
return(res)}


