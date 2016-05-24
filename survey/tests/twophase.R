library(survey)

## two-phase simple random sampling.
data(pbc, package="survival")
pbc$id<-1:nrow(pbc)
pbc$randomized<-with(pbc, !is.na(trt) & trt>-9)
(d2pbc<-twophase(id=list(~id,~id), data=pbc, subset=~I(!randomized)))
m<-svymean(~bili, d2pbc)
all.equal(as.vector(coef(m)),with(pbc, mean(bili[!randomized])))
all.equal(as.vector(SE(m)),
          with(pbc, sd(bili[!randomized])/sqrt(sum(!randomized))),
          tolerance=0.001)

## two-stage sampling as two-phase
data(mu284)
ii<-with(mu284, c(1:15, rep(1:5,n2[1:5]-3)))
mu284.1<-mu284[ii,]
mu284.1$id<-1:nrow(mu284.1)
mu284.1$sub<-rep(c(TRUE,FALSE),c(15,34-15))
dmu284<-svydesign(id=~id1+id2,fpc=~n1+n2, data=mu284)
## first phase cluster sample, second phase stratified within cluster
(d2mu284<-twophase(id=list(~id1,~id),strata=list(NULL,~id1),
                   fpc=list(~n1,NULL),data=mu284.1,subset=~sub,method="approx"))
(d22mu284<-twophase(id=list(~id1,~id),strata=list(NULL,~id1),
                   fpc=list(~n1,NULL),data=mu284.1,subset=~sub,method="full"))
summary(d2mu284)
t1<-svytotal(~y1, dmu284)
t2<-svytotal(~y1, d2mu284)
t22<-svytotal(~y1,d22mu284)
m1<-svymean(~y1, dmu284)
m2<-svymean(~y1, d2mu284)
m22<-svymean(~y1, d22mu284)
all.equal(coef(t1),coef(t2))
all.equal(coef(t1),coef(t22))
all.equal(coef(m1),coef(m2))
all.equal(coef(m1),coef(m22))
all.equal(as.vector(SE(m1)),as.vector(SE(m2)))
all.equal(as.vector(SE(m1)),as.vector(SE(m22)))
all.equal(as.vector(SE(t1)),as.vector(SE(t2)))
all.equal(as.vector(SE(t1)),as.vector(SE(t22)))

## case-cohort design
##this example requires R 2.3.1 or later for cch and data.
library("survival")
data(nwtco, package="survival")
## unstratified, equivalent to Lin & Ying (1993)
print(dcchs<-twophase(id=list(~seqno,~seqno), strata=list(NULL,~rel),
                 subset=~I(in.subcohort | rel), data=nwtco))
cch1<-svycoxph(Surv(edrel,rel)~factor(stage)+factor(histol)+I(age/12),
               design=dcchs)
dcchs2<-twophase(id=list(~seqno,~seqno), strata=list(NULL,~rel),
                 subset=~I(in.subcohort | rel), data=nwtco,method="approx")
cch1.2<-svycoxph(Surv(edrel,rel)~factor(stage)+factor(histol)+I(age/12),
               design=dcchs)
all.equal(coef(cch1),coef(cch1.2))
all.equal(SE(cch1),SE(cch1.2))
## Using survival::cch 
subcoh <- nwtco$in.subcohort
selccoh <- with(nwtco, rel==1|subcoh==1)
ccoh.data <- nwtco[selccoh,]
ccoh.data$subcohort <- subcoh[selccoh]
cch2<-cch(Surv(edrel, rel) ~ factor(stage) + factor(histol) + I(age/12),
          data =ccoh.data, subcoh = ~subcohort, id=~seqno,
          cohort.size=4028, method="LinYing", robust=TRUE)

print(all.equal(as.vector(coef(cch1)),as.vector(coef(cch2))))
## cch has smaller variances by a factor of 1.0005 because
## there is a (n/(n-1)) in the survey phase1 varianace
print(all.equal(as.vector(SE(cch1)), as.vector(SE(cch2)),tolerance=0.0006))


## bug report from Takahiro Tsuchiya for version 3.4
## We used to not match Sarndal exactly, because our old phase-one
## estimator had a small bias for finite populations
rei<-read.table(tmp<-textConnection(
"  id   N n.a h n.ah n.h   sub  y
1   1 300  20 1   12   5  TRUE  1
2   2 300  20 1   12   5  TRUE  2
3   3 300  20 1   12   5  TRUE  3
4   4 300  20 1   12   5  TRUE  4
5   5 300  20 1   12   5  TRUE  5
6   6 300  20 1   12   5 FALSE NA
7   7 300  20 1   12   5 FALSE NA
8   8 300  20 1   12   5 FALSE NA
9   9 300  20 1   12   5 FALSE NA
10 10 300  20 1   12   5 FALSE NA
11 11 300  20 1   12   5 FALSE NA
12 12 300  20 1   12   5 FALSE NA
13 13 300  20 2    8   3  TRUE  6
14 14 300  20 2    8   3  TRUE  7
15 15 300  20 2    8   3  TRUE  8
16 16 300  20 2    8   3 FALSE NA
17 17 300  20 2    8   3 FALSE NA
18 18 300  20 2    8   3 FALSE NA
19 19 300  20 2    8   3 FALSE NA
20 20 300  20 2    8   3 FALSE NA
"), header=TRUE)
close(tmp)

des.rei <- twophase(id=list(~id,~id), strata=list(NULL,~h),
                    fpc=list(~N,NULL), subset=~sub, data=rei, method="approx")
tot<- svytotal(~y, des.rei)
des.rei2 <- twophase(id=list(~id,~id), strata=list(NULL,~h),
                    fpc=list(~N,NULL), subset=~sub, data=rei)
tot2<- svytotal(~y, des.rei2)

## based on Sarndal et al (9.4.14)
rei$w.ah <- rei$n.ah / rei$n.a
a.rei <- aggregate(rei, by=list(rei$h), mean, na.rm=TRUE)
a.rei$S.ysh <- tapply(rei$y, rei$h, var, na.rm=TRUE)
a.rei$y.u <- sum(a.rei$w.ah * a.rei$y)
V <- with(a.rei, sum(N * (N-1) * ((n.ah-1)/(n.a-1) - (n.h-1)/(N-1)) * w.ah * S.ysh / n.h))
V <- V + with(a.rei, sum(N * (N-n.a) * w.ah * (y - y.u)^2 / (n.a-1)))

a.rei$f.h<-with(a.rei, n.h/n.ah)
Vphase2<-with(a.rei, sum(N*N*w.ah^2* ((1-f.h)/n.h) *S.ysh))

a.rei$f<-with(a.rei, n.a/N)
a.rei$delta.h<-with(a.rei, (1/n.h)*(n.a-n.ah)/(n.a-1))
Vphase1<-with(a.rei, sum(N*N*((1-f)/n.a)*( w.ah*(1-delta.h)*S.ysh+ ((n.a)/(n.a-1))*w.ah*(y-y.u)^2)))

V
Vphase1
Vphase2
vcov(tot)
vcov(tot2)
## phase 2 identical
all.equal(Vphase2,drop(attr(vcov(tot),"phases")$phase2))
all.equal(Vphase2,drop(attr(vcov(tot2),"phases")$phase2))
## phase 1 differs by 2.6% for old twophase estimator
Vphase1/attr(vcov(tot),"phases")$phase1
all.equal(Vphase1,as.vector(attr(vcov(tot2),"phases")$phase1))

