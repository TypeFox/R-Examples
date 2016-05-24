
## examples from Breslow & Chatterjee: Applied Statistics 1999 No. 4, p458
## data from Norman Breslow's web page.
library(survey)
load("nwts.rda")
nwtsnb<-nwts
nwtsnb$case<-nwts$case-nwtsb$case
nwtsnb$control<-nwts$control-nwtsb$control

a<-rbind(nwtsb,nwtsnb)
a$in.ccs<-rep(c(TRUE,FALSE),each=16)

b<-rbind(a,a)
b$rel<-rep(c(1,0),each=32)
b$n<-ifelse(b$rel,b$case,b$control)

index<-rep(1:64,b$n)

nwt.exp<-b[index,c(1:3,6,7)]
nwt.exp$id<-1:4088

dccs2<-twophase(id=list(~id,~id),subset=~in.ccs,
                strata=list(NULL,~interaction(instit,rel)),data=nwt.exp)

dccs8<-twophase(id=list(~id,~id),subset=~in.ccs,
                strata=list(NULL,~interaction(instit,stage,rel)),data=nwt.exp)

gccs8<-calibrate(dccs2,phase=2,formula=~interaction(instit,stage,rel))

summary(svyglm(rel~factor(stage)*factor(histol),family=quasibinomial,design=dccs2))
summary(svyglm(rel~factor(stage)*factor(histol),family=quasibinomial,design=dccs8))
summary(svyglm(rel~factor(stage)*factor(histol),family=quasibinomial,design=gccs8))

## check subsets of calibrated designs.
summary(svyglm(rel~factor(stage),
               family=quasibinomial,design=subset(dccs8,histol==1)))
summary(svyglm(rel~factor(stage),
               family=quasibinomial,design=subset(gccs8,histol==1)))

