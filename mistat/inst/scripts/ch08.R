###################################################
### Chap08Start
###################################################
library(mistat)
library(qcc)

qcc.options(bg.margin = "white",
            violating.runs = list(pch=15, col=1),
            beyond.limits = list(pch=15, col=1))



###################################################
### PlotRunChart50PistonCycleTimes
###################################################
Ps <- pistonSimulation(seed=123)

plot(Ps$seconds, 
     type="b", 
     ylab="sec.")


###################################################
### PlotHist50PistonCycelTimes
###################################################
hist(Ps$seconds, 
     main = "", 
     xlab="sec.")


###################################################
### PlotqqNorm50PistonCycleTimes
###################################################
library(car)

qqPlot(Ps$seconds, 
       distribution="norm", 
       col.lines=1, 
       ylab = "sec.", 
       envelope=FALSE)


###################################################
### Turb2Example
###################################################
Ps <- pistonSimulation(seed=123) 





Ps <- simulationGroup(Ps, 5)

head(Ps, 3)                       

aggregate(x=Ps["seconds"],        
          by=Ps["group"], 
          FUN=mean)


###################################################
### PlotXbarChart50PistonCycleTimes
###################################################
Ps <- pistonSimulation(seed=123, 
                       each=5*20)

Ps <- simulationGroup(Ps, 5)

CycleTime <- qcc.groups(Ps$seconds, 
                        Ps$group)

PsXbar <- invisible(qcc(CycleTime, 
                        type="xbar"))




###################################################
### PlotSChart50PistonCycleTimes
###################################################
PsS <- invisible(qcc(CycleTime, 
                     type="S"))

rm(CycleTime, Ps)


###################################################
### PlotXbarChart50PistonCycleTimesTrendAmbientTemperature
###################################################
PsTtm <- pistonSimulation(m = rep(60, 105),
                          s = rep(0.02, 105),
                          v0 = rep(0.01, 105),
                          k = rep(5000, 105),
                          p0 = rep(110000, 105),
                          t = c(rep(296,45), 
                                296*1.1^(1:60)),
                          t0 = rep(360, 105),
                          each = 1, 
                          seed = 123, 
                          check = FALSE)

PsTtm <- simulationGroup(PsTtm, 5)

CycleTime <- qcc.groups(PsTtm$seconds, 
                        PsTtm$group)
invisible(qcc(CycleTime, 
              type="xbar", 
              center=PsXbar$center, 
              limits=PsXbar$limits))


###################################################
### PlotSChart50PistonCycleTimesTrendAmbientTemperature
###################################################
invisible(qcc(CycleTime, 
              type="S", 
              center=PsS$center, 
              limits=PsS$limits))

rm(CycleTime, PsTtm)


###################################################
### PlotXbarChart50PistonCycleTimesTrendSpringCoefficient
###################################################
PsTk <- pistonSimulation(m = rep(60, 105),
                         s = rep(0.02, 105),
                         v0 = rep(0.01, 105),
                         k = c(rep(5000, 45), 
                               5000*0.985^(1:60)),
                         p0 = rep(110000, 105),
                         t = rep(296,105),
                         t0 = rep(360, 105),
                         each = 1, 
                         seed = 123, 
                         check = FALSE)

PsTk <- simulationGroup(PsTk, 5)

CycleTime <- qcc.groups(PsTk$seconds, 
                        PsTk$group)
invisible(qcc(CycleTime, 
              type="xbar", 
              center=PsXbar$center, 
              limits=PsXbar$limits))


###################################################
### PlotSChart50PistonCycleTimesSpringCoefficient
###################################################
invisible(qcc(CycleTime, 
              type="S", 
              center=PsS$center, 
              limits=PsS$limits))

rm(CycleTime, PsXbar, PsS, PsTk)


###################################################
### PlotXbarPatternsToDetectSpecialCauses00
###################################################
set.seed(123)

X <- rnorm(100, 10, 1)

G <- rep(1:20, each=5)

X[51:55] <- X[51:55]*1.3

Process <- qcc.groups(X, G)

layout(matrix(1:2, 2, byrow=TRUE))

invisible(qcc(Process, 
              type="xbar", 
              center=10, 
              std.dev=2.236071))

set.seed(123)

X <- rnorm(100, 10, 1)

X[21:55] <- X[21:55]*1.05

Process <- qcc.groups(X, G)

invisible(qcc(Process, 
              type="xbar", 
              center=10, 
              std.dev=2.236071))

layout(1)


###################################################
### PlotXbarPatternsToDetectSpecialCauses01
###################################################
set.seed(123)

X <- rnorm(100, 10, 1)

X[21:50] <- seq(12.5, 9.5, length.out=30)

Process <- qcc.groups(X, G)

layout(matrix(1:2, 2, byrow=TRUE))

invisible(qcc(Process, 
              type="xbar", 
              center=10, 
              std.dev=2.236071))

set.seed(123)

X <- rnorm(100, 10, 1)

X[c(21:25, 31:35)] <- 12.5

Process <- qcc.groups(X, G)

invisible(qcc(Process, 
              type="xbar", 
              center=10, 
              std.dev=2.236071))

points(c(5,7), c(12.5, 12.5), pch=15)

abline(h=10+2*2.236071/sqrt(5), lty="dotdash")

layout(1)

rm(Process, X, G)


###################################################
### Capability
###################################################
Ps <- pistonSimulation(seed=123)

Ps <- simulationGroup(Ps, 5)

CycleTime <- qcc.groups(data=Ps$seconds, 
                        sample=Ps$group)

PsXbar <- qcc(CycleTime, 
              type="xbar", 
              nsigmas=3, 
              plot=FALSE)

process.capability(PsXbar, 
                   spec.limits=c(0.1, 0.5))


###################################################
### PlotParetoChartSoftwareErrors
###################################################
data(PBX)

invisible(
  pareto.chart(PBX, 
               col=gray.colors(
                 n=length(PBX)), 
               main ="PBX software errors"))


###################################################
### PlotXbarChartShewartControlChart
###################################################
data(GASTURBINE)

invisible(qcc(GASTURBINE))

abline(h=0.4508481 + 0.04601819*2/sqrt(5), lty="dotdash")

abline(h=0.4508481 - 0.04601819*2/sqrt(5), lty="dotdash")


###################################################
### PlotPChartJanuaryData
###################################################
data(JANDEFECT)

invisible(qcc(JANDEFECT, 
              type="p", 
              sizes=100, 
              center=0.048,
              std.dev=sqrt(0.048*(1-0.048))))


###################################################
### PlotXbarChartContactLengths
###################################################
data(CONTACTLEN)

invisible(qcc(CONTACTLEN, type="xbar"))


###################################################
### PlotSChartContactLengths
###################################################
invisible(qcc(CONTACTLEN, type="S"))


###################################################
### PlotRChartContactLengths
###################################################
invisible(qcc(CONTACTLEN, type="R"))


###################################################
### Chap08End
###################################################
rm(CONTACTLEN, CycleTime, GASTURBINE, Ps, 
   JANDEFECT, PBX, PsXbar)
detach(package:qcc)
detach(package:car)
detach(package:mistat)
