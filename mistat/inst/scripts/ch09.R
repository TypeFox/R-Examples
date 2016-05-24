###################################################
### Chap09Start
###################################################
library(mistat)
library(qcc)

qcc.options(bg.margin = "white",
            violating.runs = list(pch=15, col=1),
            beyond.limits = list(pch=15, col=1))



###################################################
### PlotRandomNormalSequenceLevel10
###################################################
data(RNORM10)

plot(RNORM10, 
     type="b", 
     ylab="Sample")

abline(h=10)


###################################################
### RunsTest
###################################################
data(RNORM10)

X <- ifelse(RNORM10 <= 10, 
            yes="l", 
            no="u")

X <- as.factor(X)

library(tseries)

runs.test(X, 
          alternative="less")

rm(X)


###################################################
### PlotModifiedXbarChart50PistonCycleTimes
###################################################
Ps <- pistonSimulation(seed=123, 
                       each=5*20)

Ps <- simulationGroup(Ps, 5)

CycleTime <- qcc.groups(Ps$seconds, 
                        Ps$group)

PsXbar <- invisible(
  qcc(CycleTime, 
      type="xbar"))

St <- PsXbar$std.dev / sqrt(PsXbar$sizes[1])

abline(h=PsXbar$center + 2 * St, 
       lty="dotdash")

abline(h=PsXbar$center - 2 * St, 
       lty="dotdash")

rm(Ps, CycleTime, PsXbar, St)


###################################################
### PlotOCcurvePChart
###################################################
data(JANDEFECT)

NumDefPChart <- invisible(
  qcc(JANDEFECT, 
      type="p", 
      sizes=20, 
      center=0.048,
      std.dev=sqrt(0.048*(1-0.048)), 
      plot=FALSE))

oc.curves(NumDefPChart)
rm(NumDefPChart)


###################################################
### PlotCumulativeControlChart
###################################################
set.seed(123)

X <- c(rnorm(20, 10), rnorm(20, 13))

Cusum <- invisible(
  cusum(X, 
        center=10, 
        plot=FALSE))

plot(Cusum$pos, 
     xlab="Group", 
     ylab="Cumulative Sum", 
     type="b")

rm(Cusum, X)


###################################################
### PlotCumulativeControlChartIPL
###################################################
data(IPL)

Cusum <- cusum(IPL, 
               center=1.07, 
               std.dev=1, 
               se.shift=0, 
               sizes=1, 
               decision.interval=4.16, 
               plot=FALSE)

plot(Cusum$pos[1:30], 
     xlab="Group", 
     ylab="Cumulative Sum", 
     type="b")

abline(h=Cusum$decision.interval)

rm(Cusum)


###################################################
### PlotYearlyCoalMineDisasters
###################################################
data(COAL)                           

Bp <- barplot(COAL)                  

axis(side=1,                         
     labels=seq(
       from=1850, 
       to=1960, 
       by=10),  
     at=Bp[rep(c(TRUE,               
                 rep(FALSE, 9)),     
               10)])                 

rm(Bp)


###################################################
### PlotCusumYearlyCoalMineDisasters
###################################################
Cusum <- cusum(COAL[1:50], 
               center=1.82, 
               std.dev=1, 
               se.shift=0, 
               decision.interval=4.19, 
               plot=FALSE)

plot(Cusum$neg, 
     xlab="Group", 
     ylab="Cumulative Sum", 
     type="b")

abline(h=-Cusum$decision.interval)

rm(Cusum)


###################################################
### PlotCusumThicknessDifference
###################################################
data(THICKDIFF)

invisible(
  cusum(THICKDIFF, 
        center=0, 
        std.dev=1, 
        se.shift=6, 
        decision.interval=9, 
        plot=TRUE))


###################################################
### TableArlNorm
###################################################
cusumArl(mean= 0.0, 
         N=100,  
         limit=5000,
         seed=123)

cusumArl(mean= 0.5, 
         N=100, 
         limit=5000,
         seed=123)

cusumArl(mean= 1.0, 
         N=100,  
         limit=5000,
         seed=123)

cusumArl(mean= 1.5, 
         N=100, 
         limit=5000,
         seed=123)


###################################################
### PlotCusumRunLength
###################################################
Rl <- cusumArl(mean= 10, 
               sd= 5, 
               N=300,  
               limit=7000, 
               seed=123,
               kp=12, 
               km=8, 
               hp=29, 
               hm= -29, 
               printSummary=FALSE)

hist(Rl$rls, 
     main="", 
     xlab="Run Length", 
     breaks=18)

rm(Rl)


###################################################
### TableArlBinom
###################################################
cusumArl(size=100, 
         prob=0.05, 
         kp=5.95, 
         km=3.92, 
         hp=12.87, 
         hm=-8.66, 
         randFunc=rbinom,  
         N=100, 
         limit=2000,
         seed=123)


###################################################
### TableArlPois
###################################################
cusumArl(lambda=10, 
         kp=12.33, 
         km=8.41, 
         hp=11.36, 
         hm=-12.91, 
         randFunc=rpois, 
         N=100, 
         limit=2000,
         seed=123)


###################################################
### TablePfaCedNorm
###################################################
cusumPfaCedNorm(mean1=0.5, 
                tau=100,  
                N=100, 
                limit=1000,
                seed=123)

cusumPfaCedNorm(mean1=1.0, 
                tau=100, 
                N=100, 
                limit=1000,
                seed=123)

cusumPfaCedNorm(mean1=1.5, 
                tau=100, 
                N=100, 
                limit=1000,
                seed=123)


###################################################
### TableShroPfaCed
###################################################
shroArlPfaCedNorm(mean0=10, 
                  sd=3, 
                  n=5, 
                  delta=0.5, 
                  tau=10, 
                  w=99, 
                  seed=123)

shroArlPfaCedNorm(mean0=10, 
                  sd=3, 
                  n=5, 
                  delta=1.0, 
                  tau=10, 
                  w=99, 
                  seed=123)

shroArlPfaCedNorm(mean0=10, 
                  sd=3, 
                  n=5, 
                  delta=1.5, 
                  tau=10, 
                  w=99, 
                  seed=123)

shroArlPfaCedNorm(mean0=10, 
                  sd=3, 
                  n=5, 
                  delta=2.0, 
                  tau=10, 
                  w=99, 
                  seed=123)


###################################################
### TableShroArl
###################################################
shroArlPfaCedNorm(mean0=10, 
                  sd=3, 
                  n=5, 
                  delta=2.0, 
                  w=19, 
                  seed=123)

shroArlPfaCedNorm(mean0=10, 
                  sd=3, 
                  n=5, 
                  delta=2.0, 
                  w=50, 
                  seed=123)

shroArlPfaCedNorm(mean0=10, 
                  sd=3, 
                  n=5, 
                  delta=2.0, 
                  w=99, 
                  seed=123)


###################################################
### PlotShroRunLength
###################################################
Rl <- shroArlPfaCedNorm(mean0=10, 
                        sd=3, 
                        n=5, 
                        delta=2.0, 
                        w=99, 
                        seed=123, 
                        returnDetails=TRUE)

boxplot(Rl$rls, 
        main="", 
        ylab="Run Length")

rm(Rl)


###################################################
### PlotEwmaControlChart
###################################################
set.seed(123)

Length <- qcc.groups(
  c(rnorm(11*5, 10, 3), 
    rnorm(9*5, 14)), 
  rep(1:20, each=5))

invisible(
  ewma(Length, 
       center=10, 
       std.dev=3, 
       lambda=0.2, 
       nsigmas=3, 
       plot=TRUE))

rm(Length)


###################################################
### PlotDji1935
###################################################
data(DOJO1935)             

plot(DOJO1935,          
     type="b", 
     ylab="Dow Jones")


###################################################
### PlotFilmSpeed
###################################################
data(FILMSP)                           

Speed <- qcc.groups(FILMSP[1:215],     
                    rep(1:43, each=5)) 

invisible(                             
  ewma(Speed,                   
       center=105, 
       std.dev=6.53, 
       nsigmas=2))

rm(Speed)


###################################################
### Chap09End
###################################################
rm(COAL, DOJO1935, FILMSP, IPL, JANDEFECT, RNORM10, THICKDIFF)
detach(package:qcc)
detach(package:tseries)
detach(package:mistat)
