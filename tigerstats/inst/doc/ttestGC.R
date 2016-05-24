## ----include=FALSE-------------------------------------------------------
require(knitr)
opts_chunk$set(tidy=FALSE,fig.width=5.5,fig.height=3.5)

## ----message=FALSE-------------------------------------------------------
require(tigerstats)

## ------------------------------------------------------------------------
ttestGC(~fastest,data=m111survey)

## ------------------------------------------------------------------------
ttestGC(~fastest,data=m111survey,conf.level=0.90)

## ------------------------------------------------------------------------
ttestGC(~fastest,data=m111survey,mu=100,alternative="greater")

## ------------------------------------------------------------------------
ttestGC(~fastest,data=m111survey,mu=100,alternative="less")

## ------------------------------------------------------------------------
ttestGC(~fastest,data=m111survey,mu=100,
        alternative="two.sided")

## ----eval=FALSE----------------------------------------------------------
#  ttestGC(~fastest,data=m111survey,mu=100)

## ----ttesttutonemean2sided-----------------------------------------------
ttestGC(~fastest,data=m111survey,mu=100,
        alternative="two.sided",graph=TRUE)

## ------------------------------------------------------------------------
ttestGC(~fastest,data=m111survey,mu=100,
        alternative="two.sided",verbose=FALSE)

## ------------------------------------------------------------------------
ttestGC(mean=30,sd=4,n=40)

## ------------------------------------------------------------------------
ttestGC(mean=30,sd=4,n=40,
        mu=32,alternative="less")

## ------------------------------------------------------------------------
ttestGC(fastest~sex,data=m111survey,
        conf.level=0.85)

## ------------------------------------------------------------------------
ttestGC(fastest~sex,data=m111survey,
      mu=0)

## ------------------------------------------------------------------------
ttestGC(fastest~sex,data=m111survey,
      mu=0,first="male")

## ------------------------------------------------------------------------
ttestGC(mean=c(32,30),sd=c(4.2,5.1),n=c(33,42))

## ----ttesttut2meanssummary-----------------------------------------------
ttestGC(mean=c(32,30),sd=c(4.2,5.1),n=c(33,42),
        mu=0,alternative="greater",
        conf.level=0.90,graph=TRUE)

## ------------------------------------------------------------------------
ttestGC(~ideal_ht - height,data=m111survey)

## ------------------------------------------------------------------------
ttestGC(~ideal_ht - height,data=m111survey,
        mu=0,alternative="greater")

## ------------------------------------------------------------------------
ttestGC(~diff.ideal.act.,data=m111survey)

## ------------------------------------------------------------------------
ttestGC(~diff.ideal.act.,data=m111survey,
        mu=0,alternative="greater")

