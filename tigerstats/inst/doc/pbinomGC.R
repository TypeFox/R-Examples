## ----include=FALSE-------------------------------------------------------
require(knitr)
opts_chunk$set(tidy=FALSE,fig.width=5.5,fig.height=3.5)

## ----message=FALSE-------------------------------------------------------
require(tigerstats)

## ------------------------------------------------------------------------
pbinomGC(bound=8,region="below",size=20,prob=0.3)

## ----binomtutleq8--------------------------------------------------------
pbinomGC(bound=8,region="below",
         size=20,prob=0.3,graph=TRUE)

## ----eval=FALSE----------------------------------------------------------
#  pbinomGC(bound=10,region="above",
#           size=20,prob=0.3)

## ----binomtutgeq10-------------------------------------------------------
pbinomGC(bound=9,region="above",
         size=20,prob=0.3,graph=TRUE)

## ----binomtutless6-------------------------------------------------------
pbinomGC(bound=5,region="below",
         size=20,prob=0.3,graph=TRUE)

## ----binomtutgreater7----------------------------------------------------
pbinomGC(bound=7,region="above",
         size=20,prob=0.3,graph=TRUE)

## ----binomtutbetwee------------------------------------------------------
pbinomGC(bound=c(6,10),region="between",
         size=20,prob=0.3,graph=TRUE)

## ----binomtutequal-------------------------------------------------------
pbinomGC(bound=c(6,6),region="between",
         size=20,prob=0.3,graph=TRUE)

## ----binomtutout---------------------------------------------------------
pbinomGC(bound=c(7,9),region="outside",
         size=20,prob=0.3,graph=TRUE)

