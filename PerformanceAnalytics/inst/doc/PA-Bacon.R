### R code from vignette source 'PA-Bacon.Rnw'

###################################################
### code chunk number 1: PA-Bacon.Rnw:39-40
###################################################
library(PerformanceAnalytics)


###################################################
### code chunk number 2: PA-Bacon.Rnw:53-55
###################################################
data(portfolio_bacon)
print(MeanAbsoluteDeviation(portfolio_bacon[,1])) #expected 0.0310


###################################################
### code chunk number 3: PA-Bacon.Rnw:62-64
###################################################
data(portfolio_bacon)
print(Frequency(portfolio_bacon[,1])) #expected 12


###################################################
### code chunk number 4: PA-Bacon.Rnw:73-75
###################################################
data(managers)
SharpeRatio(managers[,1,drop=FALSE], Rf=.035/12, FUN="StdDev") 


###################################################
### code chunk number 5: PA-Bacon.Rnw:86-88
###################################################
data(portfolio_bacon)
print(MSquared(portfolio_bacon[,1], portfolio_bacon[,2])) #expected 0.1068


###################################################
### code chunk number 6: PA-Bacon.Rnw:100-103
###################################################
data(portfolio_bacon)
print(MSquaredExcess(portfolio_bacon[,1], portfolio_bacon[,2])) #expected -0.00998
print(MSquaredExcess(portfolio_bacon[,1], portfolio_bacon[,2], Method="arithmetic")) #expected -0.011


###################################################
### code chunk number 7: PA-Bacon.Rnw:116-118
###################################################
data(managers)
print(CAPM.alpha(managers[,1,drop=FALSE], managers[,8,drop=FALSE], Rf=.035/12))


###################################################
### code chunk number 8: PA-Bacon.Rnw:125-127
###################################################
data(managers)
CAPM.beta(managers[, "HAM2", drop=FALSE], managers[, "SP500 TR", drop=FALSE], Rf = managers[, "US 3m TR", drop=FALSE])


###################################################
### code chunk number 9: PA-Bacon.Rnw:138-140
###################################################
data(managers)
print(CAPM.epsilon(portfolio_bacon[,1], portfolio_bacon[,2])) #expected -0.013


###################################################
### code chunk number 10: PA-Bacon.Rnw:151-153
###################################################
data(portfolio_bacon)
print(CAPM.jensenAlpha(portfolio_bacon[,1], portfolio_bacon[,2])) #expected -0.014


###################################################
### code chunk number 11: PA-Bacon.Rnw:164-166
###################################################
data(portfolio_bacon)
print(SystematicRisk(portfolio_bacon[,1], portfolio_bacon[,2])) #expected 0.013


###################################################
### code chunk number 12: PA-Bacon.Rnw:173-175
###################################################
data(portfolio_bacon)
print(SpecificRisk(portfolio_bacon[,1], portfolio_bacon[,2])) #expected 0.0329


###################################################
### code chunk number 13: PA-Bacon.Rnw:184-186
###################################################
data(portfolio_bacon)
print(TotalRisk(portfolio_bacon[,1], portfolio_bacon[,2])) #expected 0.0134


###################################################
### code chunk number 14: PA-Bacon.Rnw:195-197
###################################################
data(managers)
print(round(TreynorRatio(managers[,1,drop=FALSE], managers[,8,drop=FALSE], Rf=.035/12),4))


###################################################
### code chunk number 15: PA-Bacon.Rnw:203-205
###################################################
data(portfolio_bacon) 
print(TreynorRatio(portfolio_bacon[,1], portfolio_bacon[,2], modified = TRUE)) #expected 1.677 


###################################################
### code chunk number 16: PA-Bacon.Rnw:216-218
###################################################
data(portfolio_bacon)
print(AppraisalRatio(portfolio_bacon[,1], portfolio_bacon[,2], method="appraisal")) #expected -0.430


###################################################
### code chunk number 17: PA-Bacon.Rnw:229-231
###################################################
data(portfolio_bacon)
print(AppraisalRatio(portfolio_bacon[,1], portfolio_bacon[,2], method="modified")) 


###################################################
### code chunk number 18: PA-Bacon.Rnw:242-244
###################################################
data(portfolio_bacon)
print(FamaBeta(portfolio_bacon[,1], portfolio_bacon[,2])) #expected 1.03


###################################################
### code chunk number 19: PA-Bacon.Rnw:255-257
###################################################
data(portfolio_bacon)
print(Selectivity(portfolio_bacon[,1], portfolio_bacon[,2])) #expected -0.0141


###################################################
### code chunk number 20: PA-Bacon.Rnw:270-272
###################################################
data(portfolio_bacon)
print(NetSelectivity(portfolio_bacon[,1], portfolio_bacon[,2])) #expected -0.017


###################################################
### code chunk number 21: PA-Bacon.Rnw:285-287
###################################################
data(managers)
TrackingError(managers[,1,drop=FALSE], managers[,8,drop=FALSE]) 


###################################################
### code chunk number 22: PA-Bacon.Rnw:298-300
###################################################
data(managers)
InformationRatio(managers[,"HAM1",drop=FALSE], managers[, "SP500 TR", drop=FALSE])


###################################################
### code chunk number 23: PA-Bacon.Rnw:313-315
###################################################
data(managers)
skewness(managers)


###################################################
### code chunk number 24: PA-Bacon.Rnw:324-326
###################################################
data(portfolio_bacon)
print(skewness(portfolio_bacon[,1], method="sample")) #expected -0.09


###################################################
### code chunk number 25: PA-Bacon.Rnw:337-339
###################################################
data(portfolio_bacon)
print(kurtosis(portfolio_bacon[,1], method="moment")) #expected 2.43


###################################################
### code chunk number 26: PA-Bacon.Rnw:348-350
###################################################
data(portfolio_bacon)
print(kurtosis(portfolio_bacon[,1], method="excess")) #expected -0.57


###################################################
### code chunk number 27: PA-Bacon.Rnw:359-361
###################################################
data(portfolio_bacon)
print(kurtosis(portfolio_bacon[,1], method="sample")) #expected 3.03


###################################################
### code chunk number 28: PA-Bacon.Rnw:370-372
###################################################
data(portfolio_bacon)
print(kurtosis(portfolio_bacon[,1], method="sample_excess")) #expected -0.41


###################################################
### code chunk number 29: PA-Bacon.Rnw:385-387
###################################################
data(portfolio_bacon)
print(PainIndex(portfolio_bacon[,1])) #expected 0.04


###################################################
### code chunk number 30: PA-Bacon.Rnw:394-396
###################################################
data(managers)
CalmarRatio(managers[,1,drop=FALSE])


###################################################
### code chunk number 31: PA-Bacon.Rnw:403-405
###################################################
data(managers)
SterlingRatio(managers[,1,drop=FALSE])


###################################################
### code chunk number 32: PA-Bacon.Rnw:416-418
###################################################
data(portfolio_bacon)
print(BurkeRatio(portfolio_bacon[,1])) #expected 0.74


###################################################
### code chunk number 33: PA-Bacon.Rnw:429-431
###################################################
data(portfolio_bacon)
print(BurkeRatio(portfolio_bacon[,1], modified = TRUE)) #expected 3.65


###################################################
### code chunk number 34: PA-Bacon.Rnw:442-444
###################################################
data(portfolio_bacon)
print(MartinRatio(portfolio_bacon[,1])) #expected 1.70


###################################################
### code chunk number 35: PA-Bacon.Rnw:455-457
###################################################
data(portfolio_bacon)
print(PainRatio(portfolio_bacon[,1])) #expected 2.66


###################################################
### code chunk number 36: PA-Bacon.Rnw:474-478
###################################################
data(portfolio_bacon)
MAR = 0.5
DownsideDeviation(portfolio_bacon[,1], MAR) #expected 0.493
DownsidePotential(portfolio_bacon[,1], MAR) #expected 0.491


###################################################
### code chunk number 37: PA-Bacon.Rnw:494-499
###################################################
data(portfolio_bacon)
MAR = 0.005
print(UpsideRisk(portfolio_bacon[,1], MAR, stat="risk")) #expected 0.02937
print(UpsideRisk(portfolio_bacon[,1], MAR, stat="variance")) #expected 0.08628
print(UpsideRisk(portfolio_bacon[,1], MAR, stat="potential")) #expected 0.01771


###################################################
### code chunk number 38: PA-Bacon.Rnw:510-513
###################################################
data(portfolio_bacon)
MAR = 0.005
print(DownsideFrequency(portfolio_bacon[,1], MAR)) #expected 0.458


###################################################
### code chunk number 39: PA-Bacon.Rnw:524-526
###################################################
data(portfolio_bacon)
print(BernardoLedoitRatio(portfolio_bacon[,1])) #expected 1.78


###################################################
### code chunk number 40: PA-Bacon.Rnw:539-541
###################################################
data(portfolio_bacon)
print(DRatio(portfolio_bacon[,1])) #expected 0.401


###################################################
### code chunk number 41: PA-Bacon.Rnw:554-557
###################################################
data(portfolio_bacon)
MAR = 0.005
print(OmegaSharpeRatio(portfolio_bacon[,1], MAR)) #expected 0.29


###################################################
### code chunk number 42: PA-Bacon.Rnw:568-570
###################################################
data(managers)
round(SortinoRatio(managers[, 1]),4)


###################################################
### code chunk number 43: PA-Bacon.Rnw:581-585
###################################################
data(portfolio_bacon)
MAR = 0.005
l = 2
print(Kappa(portfolio_bacon[,1], MAR, l)) #expected 0.157


###################################################
### code chunk number 44: PA-Bacon.Rnw:596-598
###################################################
data(edhec)
UpsidePotentialRatio(edhec[, 6], MAR=.05/12) #5 percent/yr MAR


###################################################
### code chunk number 45: PA-Bacon.Rnw:609-612
###################################################
data(portfolio_bacon)
MAR = 0.005
print(VolatilitySkewness(portfolio_bacon[,1], MAR, stat="volatility")) #expected 1.32


###################################################
### code chunk number 46: PA-Bacon.Rnw:623-626
###################################################
data(portfolio_bacon)
MAR = 0.005
print(VolatilitySkewness(portfolio_bacon[,1], MAR, stat="variability")) #expected 1.15


###################################################
### code chunk number 47: PA-Bacon.Rnw:637-639
###################################################
data(portfolio_bacon)
print(AdjustedSharpeRatio(portfolio_bacon[,1])) #expected 0.81


###################################################
### code chunk number 48: PA-Bacon.Rnw:652-654
###################################################
data(portfolio_bacon)
print(SkewnessKurtosisRatio(portfolio_bacon[,1])) #expected -0.034


###################################################
### code chunk number 49: PA-Bacon.Rnw:665-668
###################################################
data(portfolio_bacon)
MAR = 0.05
print(ProspectRatio(portfolio_bacon[,1], MAR)) #expected -0.134


###################################################
### code chunk number 50: PA-Bacon.Rnw:681-684
###################################################
data(portfolio_bacon)
MAR = 0.005
print(M2Sortino(portfolio_bacon[,1], portfolio_bacon[,2], MAR)) #expected 0.1035


###################################################
### code chunk number 51: PA-Bacon.Rnw:695-698
###################################################
data(portfolio_bacon)
MAR = 0.005
print(OmegaExcessReturn(portfolio_bacon[,1], portfolio_bacon[,2], MAR)) #expected 0.0805


###################################################
### code chunk number 52: PA-Bacon.Rnw:707-709
###################################################
data(managers)
table.Variability(managers[,1:8])


###################################################
### code chunk number 53: PA-Bacon.Rnw:717-719
###################################################
data(managers)
table.SpecificRisk(managers[,1:8], managers[,8])


###################################################
### code chunk number 54: PA-Bacon.Rnw:726-728
###################################################
data(managers)
table.InformationRatio(managers[,1:8], managers[,8])


###################################################
### code chunk number 55: PA-Bacon.Rnw:735-737
###################################################
data(managers)
table.Distributions(managers[,1:8])


###################################################
### code chunk number 56: PA-Bacon.Rnw:744-746
###################################################
data(managers)
table.DrawdownsRatio(managers[,1:8])


###################################################
### code chunk number 57: PA-Bacon.Rnw:753-755
###################################################
data(managers)
table.DownsideRiskRatio(managers[,1:8])


###################################################
### code chunk number 58: PA-Bacon.Rnw:762-764
###################################################
data(managers)
table.AnnualizedReturns(managers[,1:8])


