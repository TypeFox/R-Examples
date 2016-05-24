##
## Domain means can be written as ratio estimators or as regression coefficients
##
## This code checks that subsetting the design object gives the same results as
## these approaches.
##


library(survey)
data(fpc)
dfpc<-svydesign(id=~psuid,strat=~stratid,weight=~weight,data=fpc,nest=TRUE)
dsub<-subset(dfpc,x>4)
(m1<-svymean(~x,design=dsub))

## These should give the same domain estimates and standard errors
(m2<-svyby(~x,~I(x>4),design=dfpc, svymean,keep.var=TRUE))
m3<-svyglm(x~I(x>4)+0,design=dfpc)
summary(m3)
(m4<-svyratio(~I(x*(x>4)),~as.numeric(x>4), dfpc))
stopifnot(isTRUE(all.equal(SE(m2), as.vector(SE(m3)))))
stopifnot(isTRUE(all.equal(SE(m2)[2], as.vector(SE(m4)))))

## with strata
data(api)
dstrat<-svydesign(id=~1, strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
m1<-svymean(~enroll, subset(dstrat, comp.imp=="Yes"))
m2<-svyglm(enroll~comp.imp-1, dstrat)
m3<- svyratio(~I(enroll*(comp.imp=="Yes")), ~as.numeric(comp.imp=="Yes"), dstrat)
stopifnot(isTRUE(all.equal(as.vector(SE(m2)["comp.impYes"]), as.vector(SE(m1)))))
stopifnot(isTRUE( all.equal(as.vector(SE(m1)), as.vector(drop(SE(m3))))))

## with calibration
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
pop.totals<-c(`(Intercept)`=6194, stypeH=755, stypeM=1018)
(dclus1g3 <- calibrate(dclus1, ~stype+api99, c(pop.totals, api99=3914069)))

m1<-svymean(~api00, subset(dclus1g3, comp.imp=="Yes"))
m3<-svyratio(~I(api00*(comp.imp=="Yes")), ~as.numeric(comp.imp=="Yes"), dclus1g3)
m2<-svyglm(api00~comp.imp-1, dclus1g3)
stopifnot(isTRUE( all.equal(as.vector(SE(m2)["comp.impYes"]), as.vector(SE(m1)))))
stopifnot(isTRUE( all.equal(as.vector(SE(m1)), as.vector(drop(SE(m3))))))

## with raking
pop.types <- data.frame(stype=c("E","H","M"), Freq=c(4421,755,1018))
pop.schwide <- data.frame(sch.wide=c("No","Yes"), Freq=c(1072,5122))
dclus1r<-rake(dclus1, list(~stype,~sch.wide), list(pop.types, pop.schwide))
m1<-svymean(~api00, subset(dclus1r, comp.imp=="Yes"))
m2<-svyglm(api00~comp.imp-1, dclus1r)
m3<-svyratio(~I(api00*(comp.imp=="Yes")), ~as.numeric(comp.imp=="Yes"), dclus1r)
stopifnot(isTRUE( all.equal(as.vector(SE(m2)["comp.impYes"]), as.vector(SE(m1)))))
stopifnot(isTRUE( all.equal(as.vector(SE(m1)), as.vector(drop(SE(m3))))))



##
## based on bug report from Takahiro Tsuchiya for version 3.4
## 
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


des.rei2 <- twophase(id=list(~id,~id), strata=list(NULL,~h),
                    fpc=list(~N,NULL), subset=~sub, data=rei, method="full")
tot2<- svytotal(~y, subset(des.rei2, y>3))

rei$y<-rei$y*(rei$y>3)
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
vcov(tot2)

## comparing to regression
reg<-svyglm(y~I(y<4), design=des.rei2)
mn<-svymean(~y, subset(des.rei2,y>3))
all.equal(as.vector(coef(reg))[1],as.vector(coef(mn)))
all.equal(as.vector(SE(reg))[1],as.vector(SE(mn)))
vcov(mn)
vcov(reg)

