library(survey)
library(survival)

pbc2<-rbind(pbc,pbc)
pbc2$id<-rep(1:418,2)

dpbc1<-svydesign(id=~1, data=pbc)
dpbc2<-svydesign(id=~id, data=pbc2)

s1<-svykm(Surv(time,status>0)~1, subset(dpbc1, bili>6), se=TRUE)
s2<-svykm(Surv(time,status>0)~1, subset(dpbc2, bili>6), se=TRUE)

(c1<-confint(s1,(1:5)*365))
(c2<-confint(s2,(1:5)*365))
all.equal(c1, c2)

m1<-svycoxph(Surv(time,status>0)~log(bili), design=dpbc1)
m2<-svycoxph(Surv(time,status>0)~log(bili), design=dpbc2)

d<-data.frame(bili=c(5,10))
p1<-predict(m1, se=TRUE, newdata=d,type="curve")
p2<-predict(m2, se=TRUE, newdata=d,type="curve")

(pc1<-confint(p1[[1]],(1:5)*365))
(pc2<-confint(p2[[1]],(1:5)*365))
all.equal(pc1, pc2)

(q1<-quantile(p1[[2]]))
(q2<-quantile(p2[[2]]))
all.equal(q1,q2)
