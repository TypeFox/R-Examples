options(warn=2)
library(sensitivityPStrat)

data(vaccine.trial)

vaccine.trial.withNA <- vaccine.trial

set.seed(12345)
for(i in seq_len(20)) 
  vaccine.trial.withNA[sample(nrow(vaccine.trial), size=1, replace=TRUE),
                          sample(ncol(vaccine.trial), size=1)] <- NA

set.seed(12345)
ansJR<-with(vaccine.trial,
          sensitivityJR(z=treatment,s=hiv.outcome,y=logVL,
                    beta0=c(-1,-.75,-.5,-.25,0,.25,.5,.75,1),
                    beta1=c(-1,-.75,-.5,-.25,0,.25,.5,.75,1),
                    phi=c(0.95,1), selection="infected",
                    groupings=c("placebo","vaccine"),
                    N.boot=50)
         )
ansJR

stopifnot(is.list(ansJR))
stopifnot(inherits(ansJR,"sensitivity"))
stopifnot(inherits(ansJR,"sensitivity.0d"))
stopifnot(inherits(ansJR,"sensitivity.2.0d"))

set.seed(12345)
ansJR<-with(vaccine.trial.withNA,
          sensitivityJR(z=treatment,s=hiv.outcome,y=logVL,
                    beta0=c(-1,-.75,-.5,-.25,0,.25,.5,.75,1),
                    beta1=c(-1,-.75,-.5,-.25,0,.25,.5,.75,1),
                    phi=c(0.95, 1), selection="infected",
                    groupings=c("placebo","vaccine"), na.rm=TRUE,
                    N.boot=50)
         )
ansJR

stopifnot(is.list(ansJR))
stopifnot(inherits(ansJR,"sensitivity"))
stopifnot(inherits(ansJR,"sensitivity.0d"))
stopifnot(inherits(ansJR,"sensitivity.2.0d"))

