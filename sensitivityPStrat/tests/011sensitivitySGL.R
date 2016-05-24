options(warn=2)
library(sensitivityPStrat)

data(vaccine.trial)

vaccine.trial$followup.yearsPreART <- runif(nrow(vaccine.trial), 0.5, 3)
vaccine.trial.withNA <- vaccine.trial

set.seed(12345)
for(i in seq_len(20)) 
  vaccine.trial.withNA[sample(nrow(vaccine.trial), size=1, replace=TRUE),
                          sample(ncol(vaccine.trial), size=1)] <- NA

set.seed(12345)
sens.time<-with(vaccine.trial,
                sensitivitySGL(z=treatment, s=hiv.outcome, y=followup.yearsART,
                          d=ARTinitiation, beta=c(.25, 0,-.25,-.5), tau=3,
                          time.points=c(2,3), selection="infected",
                          trigger="initiated ART", groupings=c("placebo","vaccine"),
                          empty.principal.stratum=c("not infected","infected"),
                          N.boot=50)
               )
stopifnot(is.list(sens.time))
stopifnot(inherits(sens.time,"sensitivity"))
stopifnot(inherits(sens.time,"sensitivity.1d"))
stopifnot(all(c("Fas0", "Fas1", "beta", "alphahat") %in% names(sens.time)))
stopifnot(is.numeric(sens.time$alphahat))
stopifnot(is.numeric(sens.time$beta))
sens.time


set.seed(12345)
sens.time<-with(vaccine.trial.withNA,
                sensitivitySGL(z=treatment, s=hiv.outcome, y=followup.yearsART,
                          d=ARTinitiation, beta=c(.25, 0,-.25,-.5), tau=3,
                          time.points=c(2,3), selection="infected",
                          trigger="initiated ART", groupings=c("placebo","vaccine"),
                          empty.principal.stratum=c("not infected","infected"),
                          na.rm=TRUE, N.boot=50)
               )
sens.time

set.seed(12345)
sens.time<-with(vaccine.trial,
                sensitivitySGL(z=treatment, s=hiv.outcome, y=followup.yearsART,
                          d=ARTinitiation, v=followup.yearsPreART,
                          beta=c(.25, 0,-.25,-.5), tau=3, followup.time=2.5,
                          time.points=c(2,3), selection="infected",
                          trigger="initiated ART", groupings=c("placebo","vaccine"),
                          empty.principal.stratum=c("not infected","infected"),
                          N.boot=50)
               )
sens.time

set.seed(12345)
sens.time<-with(vaccine.trial.withNA,
                sensitivitySGL(z=treatment, s=hiv.outcome, y=followup.yearsART,
                          d=ARTinitiation, v=followup.yearsPreART,
                          beta=c(.25, 0,-.25,-.5), tau=3, followup.time=2.5,
                          time.points=c(2,3), selection="infected",
                          trigger="initiated ART", groupings=c("placebo","vaccine"),
                          empty.principal.stratum=c("not infected","infected"),
                          N.boot=50, na.rm=TRUE)
               )
sens.time

set.seed(12345)
sens.time<-with(vaccine.trial,
                sensitivitySGL(z=treatment, s=hiv.outcome, y=followup.yearsART,
                          d=ARTinitiation, v=followup.yearsPreART,
                          beta=c(.25, 0,-.25,-.5), tau=3, followup.time=2.5,
                          time.points=c(2,3), selection="infected",
                          trigger="initiated ART", groupings=c("placebo","vaccine"),
                          empty.principal.stratum=c("not infected","infected"),
                          ci.method="")
               )
sens.time
