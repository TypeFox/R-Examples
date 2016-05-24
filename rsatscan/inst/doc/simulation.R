## ------------------------------------------------------------------------
set.seed(42)
mygeo = expand.grid(1:10,1:10)
daysbase = 30
locid = rep(1:100, times=daysbase)
basecas = rbinom(3000, 1, .1)
day = rep(1:30, each = 100)
mycas = data.frame(locid,basecas, day)

## ------------------------------------------------------------------------
head(mygeo)
head(mycas)

## ------------------------------------------------------------------------
library("rsatscan")
td = tempdir()
write.geo(mygeo, location = td, file = "mygeo", userownames=TRUE)
write.cas(mycas, location = td, file = "mycas")

## ------------------------------------------------------------------------
invisible(ss.options(reset=TRUE))
ss.options(list(CaseFile="mycas.cas", PrecisionCaseTimes=4))
ss.options(list(StartDate="1", CoordinatesType=0, TimeAggregationUnits=4))
ss.options(list(EndDate="30", CoordinatesFile="mygeo.geo", AnalysisType=4, ModelType=2)) 
ss.options(list(UseDistanceFromCenterOption="y", MaxSpatialSizeInDistanceFromCenter=3)) 
ss.options(list(NonCompactnessPenalty=0, MaxTemporalSizeInterpretation=1, MaxTemporalSize=7))
ss.options(list(ProspectiveStartDate="30", ReportGiniClusters="n", LogRunToHistoryFile="n"))

## ------------------------------------------------------------------------
write.ss.prm(td, "mybase")
# This step omitted in compliance with CRAN policies
# Please install SaTScan and run the vignette with this and following code uncommented
# SaTScan can be downloaded from www.satscan.org, free of charge
# you will also find there fully compiled versions of this vignette with results

# mybase = satscan(td, "mybase")
# mybase$col[3:10]

## ------------------------------------------------------------------------
newday = data.frame(locid = 1:100, basecas = rbinom(100,1,.1), day = 31)
newcas = rbind(mycas,newday)
write.cas(newcas, location = td, file = "mycas")

## ------------------------------------------------------------------------
ss.options(list(EndDate="31"))
write.ss.prm(td, "day1")

# day1 = satscan(td, "day1")
# day1$col[3:10]

## ------------------------------------------------------------------------
newday = data.frame(locid = 1:100, basecas = rbinom(100,1,.1), day = 32)
newday$basecas[20] =5
newcas = rbind(mycas,newday)

write.cas(newcas, location = td, file = "mycas")

ss.options(list(EndDate="32"))
write.ss.prm(td, "day2")

# day2 = satscan(td,"day2")
# day2$col[3:10]

## ------------------------------------------------------------------------
# summary(day2)
# cat(day2$main[20:31],fill=1)

## ----, echo=FALSE, include=FALSE-----------------------------------------
#clean up!
file.remove(paste0(td,"/day1.prm"))
file.remove(paste0(td,"/day2.prm"))
file.remove(paste0(td,"/mybase.prm"))
file.remove(paste0(td,"/mycas.cas"))
file.remove(paste0(td,"/mygeo.geo"))

