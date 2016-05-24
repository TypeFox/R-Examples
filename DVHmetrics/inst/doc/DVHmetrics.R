## ----setup, echo=FALSE, include=FALSE, cache=FALSE------------------
# set global chunk options
knitr::opts_chunk$set(fig.align='center', fig.show='hold')
knitr::opts_chunk$set(tidy=FALSE, message=FALSE, warning=FALSE, comment=NA)
options(replace.assign=TRUE, useFancyQuotes=FALSE, show.signif.stars=FALSE, digits=4, width=70)

## ----cIntro, eval=FALSE---------------------------------------------
#  # install DVHmetrics from the CRAN online package repository
#  install.packages("DVHmetrics")

## ----cCmdline-------------------------------------------------------
## load DVHmetrics package - required for all following tasks
library(DVHmetrics, verbose=FALSE)

## calculate a DVH metric for built-in data
getMetric(dataMZ, metric="DMEAN", structure="HEART")

## ----cWebApp, eval=FALSE--------------------------------------------
#  vignette("DVHshiny")

## ----cReadData1a, eval=FALSE----------------------------------------
#  res <- readDVH("c:/folder/dvhFile.txt", type="Eclipse")

## ----cReadData1b, echo=TRUE-----------------------------------------
print(dataMZ)

## ----cReadData1c, echo=TRUE-----------------------------------------
print(dataMZ, verbose=TRUE)

## ----cReadData2, eval=FALSE-----------------------------------------
#  res <- readDVH("c:/folder/dvhFile*.txt", type="Cadplan")

## ----cReadData4, eval=FALSE-----------------------------------------
#  res <- readDVH(type="Eclipse")      # opens interactive file picker

## ----cReadData5, eval=FALSE-----------------------------------------
#  res <- readDVH("c:/folder/*", type="Eclipse", planInfo="doseRx")

## ----cMetrics1------------------------------------------------------
getMetric(dataMZ, metric="DMEAN")

## ----cMetrics2------------------------------------------------------
getMetric(dataMZ, metric="D5cc", structure="HEART")

## ----cMetrics3------------------------------------------------------
getMetric(dataMZ, metric="D5cc", structure="HEART", sortBy="observed")

## ----cMetrics4------------------------------------------------------
getMetric(dataMZ, metric=c("D10%", "V5Gy"),
          structure=c("AMYOC", "VALVE"),
          patID="23",
          sortBy=c("metric", "observed"),
          fixed=FALSE)

## ----cMetrics5------------------------------------------------------
getMetric(dataMZ, metric=c("DMEAN", "D5cc"), structure="HEART",
          sortBy="observed", splitBy="metric")

## ----cMetrics6------------------------------------------------------
met <- getMetric(dataMZ, metric=c("DMEAN", "D5cc"),
                 structure=c("HEART", "AOVALVE"),
                 sortBy="observed",
                 splitBy=c("structure", "metric"))
met                               # print the calculated results

## ----cMetricsSave1, eval=FALSE--------------------------------------
#  saveMetric(met, file="c:/folder/metrics.txt")

## ----cMetricsSave2, eval=FALSE--------------------------------------
#  saveMetric(met, file="c:/folder/metrics.txt", dec=",")

## ----cMetricsSave3, eval=FALSE--------------------------------------
#  saveMetric(met, file="c:/folder/metrics.txt", quote=TRUE)

## ----cDmean1--------------------------------------------------------
dmean <- getDMEAN(dataMZ[[1]])
subset(dmean, select=c(doseAvg, doseMed, doseMin, doseMax))

## ----cDmean2--------------------------------------------------------
# note that different tissues should have different parameter values,
# this is just for demonstration purposes
getEUD(dataMZ[[1]], EUDa=2)

## ----cDmean3--------------------------------------------------------
# note that different tissues should have different parameter values,
# this is just for demonstration purposes
getNTCP(dataMZ[[1]], NTCPtd50=40, NTCPm=0.6, NTCPn=0.5, NTCPtype="probit")

## ----cPointWise1----------------------------------------------------
# point-wise mean and SD for structure HEART over all patients
m <- getMeanDVH(dataMZ, fun=list(M=mean, SD=sd), byPat=FALSE, structure="HEART")
head(m)

## ----cPlots1, out.width='3in'---------------------------------------
showDVH(dataMZ, byPat=TRUE)

## ----cPlots2, out.width='3in'---------------------------------------
showDVH(dataMZ, byPat=FALSE, patID=c("P123", "P234"))

## ----cPlots3, out.width='3in'---------------------------------------
# match structures containing "VALVE" and "AMYOC"
showDVH(dataMZ, cumul=FALSE, rel=FALSE,
        structure=c("VALVE", "AMYOC"), fixed=FALSE)

## ----cPlots4, out.width='3in'---------------------------------------
# just save the diagram but don't show it
dvhPlot <- showDVH(dataMZ, structure=c("HEART", "AOVALVE", "AVNODE"),
                   rel=FALSE, thresh=0.001, show=FALSE)

## ----cPlots5, out.width='3in'---------------------------------------
# add point-wise mean DVH and 1 SD/2 SD regions
showDVH(dataMZ, structure="HEART", byPat=FALSE, addMSD=TRUE)

## ----cPlotsSave, eval=FALSE-----------------------------------------
#  saveDVH(dvhPlot, file="c:/folder/dvh.pdf", width=7, height=5)

## ----cConstrDef2, eval=FALSE----------------------------------------
#  dataConstr <- readConstraints("constraints.txt", dec=".", sep="\t")

## ----cConstrDef3, echo=TRUE-----------------------------------------
dataConstr     # show defined constraints and their scope

## ----cConstrCheck1, echo=TRUE---------------------------------------
## store result in object cc to save to file later
cc <- checkConstraint(dataMZ, constr=dataConstr)
print(cc, digits=2)            # show output with 2 decimal places

## ----cConstrCheck2, eval=FALSE--------------------------------------
#  saveConstraint(cc, file="c:/folder/constrCheck.txt")

## ----cConstrShow1, out.width='3in', echo=TRUE-----------------------
## plot relative volume
showConstraint(dataMZ, constr=dataConstr, byPat=TRUE)

## ----cConstrShow2, eval=FALSE---------------------------------------
#  ## plot absolute volume - store result in sc to save to file later
#  sc <- showConstraint(dataMZ, constr=dataConstr,
#                       byPat=FALSE, rel=FALSE)

## ----cConstrShow3, eval=FALSE---------------------------------------
#  saveDVH(sc, file="c:/folder/dvhConstraint.pdf")

## ----cBED1----------------------------------------------------------
getBED(D=50, fd=2.5, ab=c(2, 3, 4))
getEQD2(D=50, fd=2.5, ab=c(2, 3, 4))
getIsoEffD(D1=70, fd1=2, fd2=3, ab=c(3.5, 10))

## ----cBED2----------------------------------------------------------
getEQD2(D=dataMZ[[c(1, 1)]], fd=2.5, ab=3)

