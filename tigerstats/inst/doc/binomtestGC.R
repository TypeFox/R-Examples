## ----include=FALSE-------------------------------------------------------
require(knitr)
opts_chunk$set(tidy=FALSE,fig.width=5.5,fig.height=3.5)

## ----message=FALSE-------------------------------------------------------
require(tigerstats)

## ------------------------------------------------------------------------
binomtestGC(~sex,data=m111survey,
            success="female")

## ------------------------------------------------------------------------
binomtestGC(~sex,data=m111survey,
            success="female",conf.level=0.90)

## ------------------------------------------------------------------------
binomtestGC(~sex,data=m111survey,
            success="male",p=0.50,
            alternative="less")

## ------------------------------------------------------------------------
binomtestGC(~sex,data=m111survey,
            success="male",p=0.50,
            alternative="greater")

## ------------------------------------------------------------------------
binomtestGC(~sex,data=m111survey,
            success="male",p=0.50)

## ----binomtesttutsexless-------------------------------------------------
binomtestGC(~sex,data=m111survey,
            success="male",p=0.50,
            alternative="less",
            graph=TRUE)

## ------------------------------------------------------------------------
binomtestGC(x=1325,n=2500,p=0.50)

## ------------------------------------------------------------------------
binomtestGC(~sex,data=m111survey,
            success="male",p=0.50,
            verbose=FALSE)

