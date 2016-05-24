options(warn=2)
library(sensitivityPStrat)

data(vaccine.trial)

vaccine.trial.withNA <- vaccine.trial

set.seed(12345)
for(i in seq_len(20)) 
  vaccine.trial.withNA[sample(nrow(vaccine.trial), size=1, replace=TRUE),
                          sample(ncol(vaccine.trial), size=1)] <- NA

ans<-with(vaccine.trial,
          sensitivityGBH(z=treatment,s=hiv.outcome,y=logVL,
                    beta=c(0,.25,.5,.75,1,1.25,1.5),
                    selection="infected",
                    groupings=c("placebo","vaccine"),
                    empty.principal.stratum=c("not infected","infected"),
                    N.boot=50)
         )
ans

stopifnot(is.list(ans))
stopifnot(inherits(ans,"sensitivity"))
stopifnot(inherits(ans,"sensitivity.0d"))

ans<-with(vaccine.trial,
          sensitivityGBH(z=treatment,s=hiv.outcome,y=logVL,
                    beta=1,
                    selection="infected",
                    groupings=c("placebo","vaccine"),
                    empty.principal.stratum=c("not infected","infected"),
                    N.boot=50)
         )
ans

stopifnot(is.list(ans))
stopifnot(inherits(ans,"sensitivity"))
stopifnot(inherits(ans,"sensitivity.0d"))

ans<-with(vaccine.trial,
          sensitivityGBH(z=treatment,s=hiv.outcome,y=logVL,
                    beta=c(-Inf, 0,.25,.5,.75,1,Inf),                      
                    selection="infected",
                    groupings=c("placebo","vaccine"),
                    empty.principal.stratum=c("not infected","infected"),
                    ci.method="bootstrap",
                    N.boot=50)
         )
ans
stopifnot(is.list(ans))
stopifnot(inherits(ans,"sensitivity"))
stopifnot(inherits(ans,"sensitivity.0d"))

ans<-with(vaccine.trial,
          sensitivityGBH(z=treatment,s=hiv.outcome,y=logVL,
                    beta=-Inf,                      
                    selection="infected",
                    groupings=c("placebo","vaccine"),
                    empty.principal.stratum=c("not infected","infected"),
                    ci.method="bootstrap",
                    N.boot=50)
         )
ans

stopifnot(is.list(ans))
stopifnot(inherits(ans,"sensitivity"))
stopifnot(inherits(ans,"sensitivity.0d"))

ans<-with(vaccine.trial,
          sensitivityGBH(z=treatment,s=hiv.outcome,y=logVL,
                    beta=c(-Inf, seq(-5,5,length=21),Inf),
                    selection="infected",
                    groupings=c("placebo","vaccine"),
                    empty.principal.stratum=c("not infected","infected"),
                    ci.method="bootstrap",
                    N.boot=50)
         )
ans
stopifnot(is.list(ans))
stopifnot(inherits(ans,"sensitivity"))
stopifnot(inherits(ans,"sensitivity.0d"))

ans<-with(vaccine.trial,
          sensitivityGBH(z=treatment,s=hiv.outcome,y=logVL,
                    beta=c(0,.25,.5,.75,1,1.25,1.5),
                    selection="infected",
                    groupings=c("placebo","vaccine"),
                    ci.method="bootstrap",
                    method=c("ACE","T1","T2"),
                    empty.principal.stratum=c("not infected","infected"),
                    N.boot=50)
         )
ans

stopifnot(is.list(ans))
stopifnot(inherits(ans,"sensitivity"))
stopifnot(inherits(ans,"sensitivity.0d"))

ans<-with(vaccine.trial,
          sensitivityGBH(z=treatment,s=hiv.outcome,y=logVL,
                    beta=1,
                    selection="infected",
                    groupings=c("placebo","vaccine"),
                    ci.method="bootstrap",
                    method=c("ACE","T1","T2"),
                    empty.principal.stratum=c("not infected","infected"),
                    N.boot=50)
         )
ans

stopifnot(is.list(ans))
stopifnot(inherits(ans,"sensitivity"))
stopifnot(inherits(ans,"sensitivity.0d"))

ans<-with(vaccine.trial.withNA,
          sensitivityGBH(z=treatment,s=hiv.outcome,y=logVL,
                    beta=c(0,.25,.5,.75,1,1.25,1.5),
                    selection="infected",
                    groupings=c("placebo","vaccine"),
                    empty.principal.stratum=c("not infected","infected"),
                    na.rm=TRUE, N.boot=50)
         )
ans
