## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ------------------------------------------------------------------------
library(etma)

## ------------------------------------------------------------------------
data(data.GST)
head(data.GST)

## ------------------------------------------------------------------------
data(data.PAH)
head(data.PAH)

## ------------------------------------------------------------------------
data(data.RAS)
head(data.RAS)

## ----results = 'hide'----------------------------------------------------
ggint.toy=ETMA(case.ACE.0,case.ACE.1,ctrl.ACE.0,ctrl.ACE.1,
                  case.AGT.0,case.AGT.1,ctrl.AGT.0,ctrl.AGT.1,
                  data=data.RAS,iterations.step1=100,iterations.step2=300,
                  start.seed=1,show.detailed.plot=FALSE,show.final.plot=FALSE)

## ------------------------------------------------------------------------
print(ggint.toy)
summary(ggint.toy)

