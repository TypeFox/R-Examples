library(survey)
library(survival)
data(nwtco)

ntwco<-subset(nwtco, !is.na(edrel))

load("nwtco-subcohort.rda")
nwtco$subcohort<-subcohort

d_BorganII <- twophase(id=list(~seqno,~seqno),
                       strata=list(NULL,~interaction(instit,rel)),
                       data=nwtco, subset=~I(rel |subcohort))

##Coefficient results same as Splus with code from
## http://faculty.washington.edu/norm/software.html
## SE slightly larger due to using sandwich variance.

svycoxph(Surv(edrel, rel)~factor(stage)+factor(histol)+I(age/12), design=d_BorganII)

##
## This gives higher standard errors. calibrate() does not recompute the
## finite population correction if a calibration variable happens to predict
## sampling perfectly. It probably should.
##
d_BorganIIps<-calibrate(twophase(id=list(~seqno,~seqno),
                           strata=list(NULL,~rel),
                           data=nwtco, subset=~I(rel |subcohort)), 
                        phase=2, formula=~interaction(instit,rel),
                        epsilon=1e-10)

svycoxph(Surv(edrel, rel)~factor(stage)+factor(histol)+I(age/12), design=d_BorganIIps)
