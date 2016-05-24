options(warn=2)
library(sensitivityPStrat)

data(vaccine.trial)

vaccine.trial.withNA <- vaccine.trial

set.seed(12345)
for(i in seq_len(20)) 
  vaccine.trial.withNA[sample(nrow(vaccine.trial), size=1, replace=TRUE),
                          sample(ncol(vaccine.trial), size=1)] <- NA

set.seed(12345)
est.bounds<-with(vaccine.trial,
                 sensitivityHHS(z=treatment, s=hiv.outcome, y=logVL,
                     selection="infected", groupings=c("placebo","vaccine"),
                     empty.principal.stratum=c("not infected","infected"),
                     ci = 0.9, ci.method="bootstrap", ci.type="lower",
                     method=c("T1"), N.boot=50)
                )
est.bounds

set.seed(12345)
est.bounds<-with(vaccine.trial,
                 sensitivityHHS(z=treatment, s=hiv.outcome, y=logVL,
                     selection="infected", groupings=c("placebo","vaccine"),
                     empty.principal.stratum=c("not infected","infected"),
                     ci = 0.9, ci.method="bootstrap", ci.type="lower",
                     method=c("ACE", "T1", "T2"), N.boot=50)
                )
est.bounds

set.seed(12345)
est.bounds<-with(vaccine.trial,
                 sensitivityHHS(z=treatment, s=hiv.outcome, y=logVL,
                     selection="infected", groupings=c("placebo","vaccine"),
                     empty.principal.stratum=c("not infected","infected"),
                     ci = c(0.95, 0.9, 0.9), ci.method="bootstrap",
                     ci.type=c("twoSided", "lower", "upper"),
                     method=c("ACE", "T1", "T2"), N.boot=50)
                )
est.bounds

set.seed(12345)
est.bounds<-with(vaccine.trial,
                 sensitivityHHS(z=treatment, s=hiv.outcome, y=logVL,
                     selection="infected", groupings=c("placebo","vaccine"),
                     empty.principal.stratum=c("not infected","infected"),
                     ci = c(0.95, 0.9), ci.method="bootstrap",
                     ci.type=c("twoSided", "lower"),
                     method=c("T1", "T2"), N.boot=50)
                )
est.bounds

set.seed(12345)
est.bounds<-with(vaccine.trial,
                 sensitivityHHS(z=treatment, s=hiv.outcome, y=logVL,
                     selection="infected", groupings=c("placebo","vaccine"),
                     empty.principal.stratum=c("not infected","infected"),
                     method=c("ACE", "T1", "T2"), N.boot=50)
                )
est.bounds


set.seed(12345)
est.bounds<-with(vaccine.trial.withNA,
                 sensitivityHHS(z=treatment, s=hiv.outcome, y=logVL,
                     selection="infected", groupings=c("placebo","vaccine"),
                     empty.principal.stratum=c("not infected","infected"),     
                     method=c("ACE", "T1", "T2"), na.rm=TRUE, N.boot=50)
                )
est.bounds
