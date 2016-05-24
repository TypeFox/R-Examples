## ----include=FALSE-------------------------------------------------------
require(knitr)
opts_chunk$set(tidy=FALSE,fig.width=5.5,fig.height=3.5)

## ----message=FALSE-------------------------------------------------------
library(tigerstats)

## ------------------------------------------------------------------------
xtabs(~sex+seat,data=m111survey)

## ----eval=FALSE----------------------------------------------------------
#  helpGC(bwplot)

## ----eval=FALSE----------------------------------------------------------
#  helpGC(ttestGC)

## ----eval=FALSE----------------------------------------------------------
#  shiny::runApp(system.file("SlowGoodness",package="tigerstats"))

