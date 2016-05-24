## ----include=FALSE-------------------------------------------------------
require(knitr)
opts_chunk$set(tidy=FALSE,fig.width=5.5,fig.height=3.5)

## ----include=FALSE-------------------------------------------------------
library(mosaic)
library(tigerstats)

## ----eval=FALSE----------------------------------------------------------
#  binomtestGC(~sex,data=m111survey,
#              p=0.50,
#              alternative="two.sided",
#              success="female")

## ----eval=FALSE----------------------------------------------------------
#  chisqtestGC(~seat,data=m111survey,p=c(1/3,1/3,1/3))

## ----eval=FALSE----------------------------------------------------------
#  proptestGC(~sex+love_first,data=m111survey,
#                  p=0,
#                  alternative="two,sided")

## ----eval=FALSE----------------------------------------------------------
#  chisqtestGC(~sex+seat,data=m111survey)

## ----eval=FALSE----------------------------------------------------------
#  ttestGC(~fastest,data=m111survey,mu=100,
#            alternative="two.sided")

## ----eval=FALSE----------------------------------------------------------
#  ttestGC(~height - ideal_ht,data=m111survey,
#            mu=0,
#            alternative="two.sided")

## ----eval=FALSE----------------------------------------------------------
#  ttestGC(fastest~sex,data=m111survey,
#            mu=0,
#            alternative="two.sided")

