## ----include=FALSE-------------------------------------------------------
require(knitr)
opts_chunk$set(tidy=FALSE,fig.width=5.5,fig.height=3.5)

## ----message=FALSE-------------------------------------------------------
require(tigerstats)

## ------------------------------------------------------------------------
proptestGC(~sex,data=m111survey,
           success="female")

## ------------------------------------------------------------------------
proptestGC(~sex,data=m111survey,
           success="female",conf.level=0.90)

## ------------------------------------------------------------------------
proptestGC(~sex,data=m111survey,p=0.50,
           success="female",alternative="greater")

## ------------------------------------------------------------------------
proptestGC(~sex,data=m111survey,p=0.50,
           success="female",alternative="less")

## ------------------------------------------------------------------------
proptestGC(~sex,data=m111survey,p=0.50,
        success="female",alternative="two.sided")

## ----eval=FALSE----------------------------------------------------------
#  proptestGC(~sex,data=m111survey,p=0.50,
#             success="female")

## ----proptesttutoneprop2sided--------------------------------------------
proptestGC(~sex,data=m111survey,p=0.50,
        success="female",graph=TRUE)

## ------------------------------------------------------------------------
proptestGC(~sex,data=m111survey,p=0.50,
        success="female",verbose=FALSE)

## ------------------------------------------------------------------------
proptestGC(x=1325,n=2500)

## ------------------------------------------------------------------------
proptestGC(x=1325,n=2500,
        p=0.50,,alternative="less")

## ------------------------------------------------------------------------
proptestGC(~sex+enough_Sleep,data=m111survey,
        success="yes",conf.level=0.85)

## ------------------------------------------------------------------------
proptestGC(~sex+enough_Sleep,data=m111survey,
      success="yes",p=0)

## ------------------------------------------------------------------------
proptestGC(~sex+enough_Sleep,data=m111survey,
      success="yes",p=0,first="male")

## ------------------------------------------------------------------------
proptestGC(x=c(1325,905),n=c(2500,1800))

## ----proptesttut2propssummary--------------------------------------------
proptestGC(x=c(1325,905),n=c(2500,1800),
           p=0,conf.level=0.90,graph=TRUE)

