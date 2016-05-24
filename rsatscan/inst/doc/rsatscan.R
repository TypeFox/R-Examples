## ------------------------------------------------------------------------
library("rsatscan")

## ------------------------------------------------------------------------
head(NYCfevercas)
head(NYCfevergeo)

## ------------------------------------------------------------------------
invisible(ss.options(reset=TRUE))

## ------------------------------------------------------------------------
ss.options(list(CaseFile="NYCfever.cas", PrecisionCaseTimes=3))
ss.options(c("StartDate=2001/11/1","EndDate=2001/11/24"))
ss.options(list(CoordinatesFile="NYCfever.geo", AnalysisType=4, ModelType=2, TimeAggregationUnits=3))
ss.options(list(UseDistanceFromCenterOption="y", MaxSpatialSizeInDistanceFromCenter=3, NonCompactnessPenalty=0))
ss.options(list(MaxTemporalSizeInterpretation=1, MaxTemporalSize=7))
ss.options(list(ProspectiveStartDate="2001/11/24", ReportGiniClusters="n", LogRunToHistoryFile="n"))

## ------------------------------------------------------------------------
head(ss.options(),3)

## ------------------------------------------------------------------------
td = tempdir()
write.ss.prm(td, "NYCfever")
write.cas(NYCfevercas, td, "NYCfever")
write.geo(NYCfevergeo, td, "NYCfever")

## ------------------------------------------------------------------------
# This step omitted in compliance with CRAN policies
# Please install SaTScan and run the vignette with this and following code uncommented
# SaTScan can be downloaded from www.satscan.org, free of charge
# you will also find there fully compiled versions of this vignette with results

## NYCfever = satscan(td, "NYCfever", sslocation="C:/Program Files (x86)/SaTScan")

## ------------------------------------------------------------------------
## summary(NYCfever)

## ------------------------------------------------------------------------
## summary.default(NYCfever)

## ------------------------------------------------------------------------
## sp::plot(NYCfever$shapeclust)

## ----, fig.keep="all", fig.show="hold"-----------------------------------
## hist(unlist(NYCfever$llr), main="Monte Carlo")

# Let's draw a line for the clusters in the observed data
## abline(v=NYCfever$col[,c("TEST_STAT")], col = "red")

## ----, echo=FALSE, include=FALSE-----------------------------------------
#clean up!
file.remove(paste0(td,"/NYCfever.prm"))
file.remove(paste0(td,"/NYCfever.cas"))
file.remove(paste0(td,"/NYCfever.geo"))

## ------------------------------------------------------------------------
write.cas(NMcas, td,"NM")
write.geo(NMgeo, td,"NM")
write.pop(NMpop, td,"NM")

## ------------------------------------------------------------------------
invisible(ss.options(reset=TRUE))
ss.options(list(CaseFile="NM.cas",StartDate="1973/1/1",EndDate="1991/12/31", 
                PopulationFile="NM.pop",
                CoordinatesFile="NM.geo", CoordinatesType=0, AnalysisType=3))
ss.options(c("NonCompactnessPenalty=0", "ReportGiniClusters=n", "LogRunToHistoryFile=n"))

write.ss.prm(td,"testnm")
## testnm = satscan(td,"testnm", sslocation="C:/Program Files (x86)/SaTScan")

## ------------------------------------------------------------------------
## summary(testnm)

## ------------------------------------------------------------------------
## head(testnm$prm,10)

## ----, echo=FALSE, include=FALSE-----------------------------------------
#clean up!
file.remove(paste0(td,"/testnm.prm"))
file.remove(paste0(td,"/NM.pop"))
file.remove(paste0(td,"/NM.cas"))
file.remove(paste0(td,"/NM.geo"))

## ------------------------------------------------------------------------
write.cas(NHumbersidecas, td, "NHumberside")
write.ctl(NHumbersidectl, td, "NHumberside")
write.geo(NHumbersidegeo, td, "NHumberside")

invisible(ss.options(reset=TRUE))
ss.options(list(CaseFile="NHumberside.cas", ControlFile="NHumberside.ctl"))
ss.options(list(PrecisionCaseTimes=0, StartDate="2001/11/1", EndDate="2001/11/24"))
ss.options(list(CoordinatesFile="NHumberside.geo", CoordinatesType=0, ModelType=1))
ss.options(list(TimeAggregationUnits = 3, NonCompactnessPenalty=0))
ss.options(list(ReportGiniClusters="n", LogRunToHistoryFile="n"))

write.ss.prm(td, "NHumberside")
## NHumberside = satscan(td, "NHumberside", sslocation="C:/Program Files (x86)/SaTScan")

## summary(NHumberside)

## ----, echo=FALSE, include=FALSE-----------------------------------------
#clean up!
file.remove(paste0(td,"/NHumberside.cas"))
file.remove(paste0(td,"/NHumberside.ctl"))
file.remove(paste0(td,"/NHumberside.geo"))
file.remove(paste0(td,"/NHumberside.prm"))

