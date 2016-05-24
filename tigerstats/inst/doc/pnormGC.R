## ----include=FALSE-------------------------------------------------------
require(knitr)
opts_chunk$set(tidy=FALSE,fig.width=5.5,fig.height=3.5)

## ----message=FALSE-------------------------------------------------------
require(tigerstats)

## ------------------------------------------------------------------------
pnormGC(bound=66,region="below",mean=70,sd=3)

## ----normtutless66-------------------------------------------------------
pnormGC(bound=66,region="below",
        mean=70,sd=3,graph=TRUE)

## ----normtutgreater69----------------------------------------------------
pnormGC(bound=69,region="above",
        mean=70,sd=3,graph=TRUE)

## ----normtutbetween------------------------------------------------------
pnormGC(bound=c(68,72),region="between",
        mean=70,sd=3,graph=TRUE)

## ----mormtutout----------------------------------------------------------
pnormGC(bound=c(66,73),region="outside",
         mean=70,sd=3,graph=TRUE)

## ----eval=FALSE----------------------------------------------------------
#  pnormGC(bound=66,region="below",
#          mean=70,sd=3)

## ----normtutequal--------------------------------------------------------
pnormGC(bound=c(66,66),region="between",
        mean=70,sd=3,graph=TRUE)

