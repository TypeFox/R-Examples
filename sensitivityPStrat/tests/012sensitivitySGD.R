options(warn=2)
library(sensitivityPStrat)

data(vaccine.trial)

vaccine.trial$followup.yearsPreART <- runif(2000, 0.5, 3)
vaccine.trial.withNA <- vaccine.trial

set.seed(12345)
for(i in seq_len(20)) 
  vaccine.trial.withNA[sample(nrow(vaccine.trial), size=1, replace=TRUE),
                          sample(ncol(vaccine.trial), size=1)] <- NA

set.seed(12345)
sens.analysis<-with(vaccine.trial,
                sensitivitySGD(z=treatment, s=hiv.outcome, y=followup.yearsART,
                          d=ARTinitiation, beta0=c(0,-.25,-.5),
                          beta1=c(0, -.25, -.5), phi=c(0.95, 0.90, 1), tau=3,
                          time.points=c(2,3), selection="infected",
                          trigger="initiated ART",
                          groupings=c("placebo","vaccine"), ci=.95,
                          ci.method="bootstrap", N.boot=50)
               )
stopifnot(is.list(sens.analysis))
stopifnot(inherits(sens.analysis,"sensitivity"))
stopifnot(inherits(sens.analysis,"sensitivity.1d"))
stopifnot(all(c("Fas0", "Fas1", "beta0", "alphahat0", "beta1", "alphahat1") %in% names(sens.analysis)))
stopifnot(is.numeric(sens.analysis$alphahat0))
stopifnot(is.numeric(sens.analysis$beta0))
stopifnot(is.numeric(sens.analysis$alphahat1))
stopifnot(is.numeric(sens.analysis$beta1))
stopifnot(with(sens.analysis,
               Fas0[1, 1](2) - Fas1[1, 1](2) == SCE[1,1,1,1]))
sens.analysis



set.seed(12345)
sens.analysis<-with(vaccine.trial.withNA,
                sensitivitySGD(z=treatment, s=hiv.outcome, y=followup.yearsART,
                          d=ARTinitiation, beta0=c(0,-.25,-.5),
                          beta1=c(0, -.25, -.5), phi=c(0.95, 0.90, 1), tau=3,
                          time.points=c(2,3), selection="infected",
                          trigger="initiated ART",
                          groupings=c("placebo","vaccine"), ci=.95, na.rm=TRUE,
                          ci.method="bootstrap", N.boot=50)
               )
sens.analysis


set.seed(12345)
sens.analysis<-with(vaccine.trial,
                sensitivitySGD(z=treatment, s=hiv.outcome, y=followup.yearsART,
                          v=followup.yearsPreART, d=ARTinitiation,
                          beta0=c(0,-.25,-.5),
                          beta1=c(0, -.25, -.5), phi=c(0.95, 0.90, 1), tau=3,
                          followup.time=2.5,
                          time.points=c(2,3), selection="infected",
                          trigger="initiated ART",
                          groupings=c("placebo","vaccine"), ci=.95,
                          ci.method="bootstrap", N.boot=50)
               )
sens.analysis


set.seed(12345)
sens.analysis<-with(vaccine.trial.withNA,
                sensitivitySGD(z=treatment, s=hiv.outcome, y=followup.yearsART,
                          v=followup.yearsPreART, d=ARTinitiation,
                          beta0=c(0,-.25,-.5),
                          beta1=c(0, -.25, -.5), phi=c(0.95, 0.90, 1), tau=3,
                          followup.time=2.5,
                          time.points=c(2,3), selection="infected",
                          trigger="initiated ART",
                          groupings=c("placebo","vaccine"), ci=.95,
                          ci.method="bootstrap", N.boot=50, na.rm=TRUE)
               )
sens.analysis


set.seed(12345)
sens.analysis<-with(vaccine.trial,
                sensitivitySGD(z=treatment, s=hiv.outcome, y=followup.yearsART,
                          d=ARTinitiation, beta0=c(0,-.25,-.5),
                          beta1=c(0, -.25, -.5), phi=c(1), tau=3,
                          time.points=c(2,3), selection="infected",
                          trigger="initiated ART",
                          groupings=c("placebo","vaccine"), ci=.95,
                          ci.method="")
               )
sens.analysis
