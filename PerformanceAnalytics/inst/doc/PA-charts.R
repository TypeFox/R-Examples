### R code from vignette source 'PA-charts.Rnw'

###################################################
### code chunk number 1: LoadLibrary
###################################################
library('PerformanceAnalytics')
#source("PerformanceAnalytics/R/charts.PerformanceSummary.R")


###################################################
### code chunk number 2: LoadData
###################################################
data(managers)
#managers=read.csv("/home/peter/dev/R/managers.csv",row.names=1)
head(managers)


###################################################
### code chunk number 3: CalcDataDimensions
###################################################
dim(managers)
managers.length = dim(managers)[1]
colnames(managers)
manager.col = 1
peers.cols = c(2,3,4,5,6)
indexes.cols = c(7,8)
Rf.col = 10
#factors.cols = NA
trailing12.rows = ((managers.length - 11):managers.length)
trailing12.rows
trailing36.rows = ((managers.length - 35):managers.length)
trailing60.rows = ((managers.length - 59):managers.length)
#assume contiguous NAs - this may not be the way to do it na.contiguous()?
frInception.rows = (length(managers[,1]) - length(managers[,1][!is.na(managers[,1])]) + 1):length(managers[,1])



###################################################
### code chunk number 4: ShowColors
###################################################
tim12equal


###################################################
### code chunk number 5: ShowPalettes
###################################################
layout(rbind(c(1,2),c(3,4)))
chart.CumReturns(managers[,c(manager.col,peers.cols)],
lwd=2,legend.loc="topleft", main="default")
chart.CumReturns(managers[,c(manager.col,peers.cols)],lwd=2,
legend.loc="topleft",colorset=redfocus, main = "redfocus")
chart.CumReturns(managers[,c(manager.col,peers.cols)],lwd=2,
legend.loc="topleft",colorset=rainbow6equal, main="rainbow6equal")
chart.CumReturns(managers[,c(manager.col,peers.cols)],lwd=2,
legend.loc="topleft",colorset=tim6equal, main = "tim6equal")


###################################################
### code chunk number 6: Graph1
###################################################
charts.PerformanceSummary(managers[,c(manager.col,indexes.cols)], colorset=rich6equal, lwd=2, ylog=TRUE)


###################################################
### code chunk number 7: CalendarReturns
###################################################
t(table.CalendarReturns(managers[,c(manager.col,indexes.cols)]))


###################################################
### code chunk number 8: MonthlyReturnStats
###################################################
table.Stats(managers[,c(manager.col,peers.cols)])


###################################################
### code chunk number 9: Graph10
###################################################
chart.Boxplot(managers[ trailing36.rows, c(manager.col, peers.cols, indexes.cols)], main = "Trailing 36-Month Returns")


###################################################
### code chunk number 10: Graph13
###################################################
layout(rbind(c(1,2),c(3,4)))
chart.Histogram(managers[,1,drop=F], main = "Plain", methods = NULL)
chart.Histogram(managers[,1,drop=F], main = "Density", breaks=40,
methods = c("add.density", "add.normal"))
chart.Histogram(managers[,1,drop=F], main = "Skew and Kurt", methods = c
("add.centered", "add.rug"))
chart.Histogram(managers[,1,drop=F], main = "Risk Measures", methods = c
("add.risk"))


###################################################
### code chunk number 11: Graph3
###################################################
chart.RiskReturnScatter(managers[trailing36.rows,1:8], Rf=.03/12, main = "Trailing 36-Month Performance", colorset=c("red", rep("black",5), "orange", "green"))


###################################################
### code chunk number 12: Graph5
###################################################
charts.RollingPerformance(managers[, c(manager.col, peers.cols, indexes.cols)], Rf=.03/12, colorset = c("red", rep("darkgray",5), "orange", "green"), lwd = 2)


###################################################
### code chunk number 13: Graph6
###################################################
chart.RelativePerformance(managers[ , manager.col, drop = FALSE], managers[ , c(peers.cols, 7)], colorset = tim8equal[-1], lwd = 2, legend.loc = "topleft")


###################################################
### code chunk number 14: Graph6a
###################################################
chart.RelativePerformance(managers[ , c(manager.col, peers.cols) ], managers[, 8, drop=F], colorset = rainbow8equal, lwd = 2, legend.loc = "topleft")


###################################################
### code chunk number 15: tableCAPM
###################################################
table.CAPM(managers[trailing36.rows, c(manager.col, peers.cols)],
managers[ trailing36.rows, 8, drop=FALSE],
Rf = managers[ trailing36.rows, Rf.col, drop=FALSE])


###################################################
### code chunk number 16: Graph8
###################################################
#source("PerformanceAnalytics/R/Return.excess.R")
charts.RollingRegression(managers[, c(manager.col, peers.cols), drop = FALSE], managers[, 8, drop = FALSE], Rf = .03/12, colorset = redfocus, lwd = 2)


###################################################
### code chunk number 17: Graph9
###################################################
chart.RollingCorrelation(managers[,c(manager.col, peers.cols)], managers[, 8, drop = FALSE], colorset = tim8equal, lwd = 2, main = "12-Month Rolling Correlation")


###################################################
### code chunk number 18: tableCorr
###################################################
table.Correlation(managers[, c(manager.col, peers.cols)], managers[, 8,
drop = F], legend.loc = "lowerleft")


###################################################
### code chunk number 19: tableDownside
###################################################
table.DownsideRisk(managers[,1:6],Rf=.03/12)


###################################################
### code chunk number 20: tableDrawdowns
###################################################
table.Drawdowns(managers[,1,drop=F])


