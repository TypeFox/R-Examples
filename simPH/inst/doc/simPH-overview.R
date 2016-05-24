
## ----include=FALSE-------------------------------------------------------
#### Load Packages ####
library(knitr)
library(simPH)
library(survival)
library(ggplot2)
library(gridExtra)

#### Load data ####
# Load Carpenter (2002) data 
## The data is included with simPH
data("CarpenterFdaData")

# Load Golub & Steunenberg (2007) data
## The data is included with simPH
data("GolubEUPData")

options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)

##### Set Chunk Options ####
opts_chunk$set(fig.align='center', dev='png', prompt=TRUE, highlight=FALSE, background="white")


## ----LinearModel1_1, tidy=FALSE, echo=TRUE, message=FALSE, warning=FALSE----
library("survival") 
library("simPH")

hmohiv <- read.table(
            "http://www.ats.ucla.edu/stat/r/examples/asa/hmohiv.csv", 
            sep = ",", header = TRUE)


## ----LinearModel1_2, tidy=FALSE, echo=TRUE, message=FALSE, warning=FALSE----
hmohiv$AgeMed <- hmohiv$age - 35 


## ----LinearModel1.3, tidy=FALSE, echo=TRUE, message=FALSE, warning=FALSE----
M1 <- coxph(Surv(time, censor) ~ AgeMed + drug,  
            method = "breslow", data = hmohiv)


## ----LinearModel1_4, tidy=FALSE, echo=TRUE, message=FALSE, warning=FALSE----
Sim1 <- coxsimLinear(M1, b = "AgeMed", Xj = seq(-15, 19, by = 1),
                    nsim = 100)


## ----LinearModel1_5, eval=FALSE------------------------------------------
## simGG(Sim1)


## ----LinearModel1_6, eval=FALSE, tidy=FALSE, echo=TRUE, message=FALSE, warning=FALSE----
## simGG(Sim1, xlab = "\nYears of Age from the Sample Median (35)",
##       ylab = "Relative Hazard with Comparison\n to a 35 Year Old\n",
##       alpha = 0.05, type = "lines")


## ----LinearModel1_7, tidy=FALSE, echo=FALSE, message=FALSE, warning=FALSE, fig.height=4, out.width='0.95\\linewidth', dev='pdf'----
# Plot results with simGG default
Plot1_1 <- simGG(Sim1)

# Plot results
Plot1_2 <-simGG(Sim1, xlab = "\nYears of Age from the\n Sample Median (35)",
                ylab = "Relative Hazard with Comparison\n to a 35 Year Old\n",
                alpha = 0.05, type = "lines")

# Combine plots
grid.arrange(Plot1_1, Plot1_2, ncol = 2)


## ----LinearModel1_8, tidy=FALSE, message=FALSE---------------------------
Sim2 <- coxsimLinear(M1, b = "drug", Xj = 0:1, nsim = 100)


## ----LinearModel1_9, tidy=FALSE, echo=FALSE, message=FALSE, warning=FALSE, fig.width=7, fig.height=4, out.width='0.6\\linewidth', dev='pdf'----
simGG(Sim2, psize = 3, xlab = "",
      ylab = "Relative Hazard\n",
      type = 'points', method = 'lm') + 
    scale_x_continuous(breaks = c(0, 1), 
                       labels = c('\nNo Drug Use', '\nDrug Use'))


## ----TVCModel1_1, tidy=FALSE, echo=TRUE, message=FALSE, warning=FALSE----
data("GolubEUPData")


## ----TVCModel1_1_1, tidy=FALSE, echo=TRUE, message=FALSE, warning=FALSE----
head(GolubEUPData[, 2:5])


## ----TVCModel1_1_2, tidy=FALSE, echo=TRUE, message=FALSE, warning=FALSE----
GolubEUPData <- SurvExpand(GolubEUPData, GroupVar = 'caseno',
                      Time = 'begin', Time2 = 'end', event = 'event') 

head(GolubEUPData[, 1:4])


## ----TVCModel1_2, tidy=FALSE, echo=TRUE, message=FALSE, warning=FALSE----
BaseVars <- c('qmv', 'backlog', 'coop', 'codec', 'qmvpostsea', 'thatcher')

GolubEUPData <- tvc(GolubEUPData, b = BaseVars, tvar = 'end', tfun = 'log')

names(GolubEUPData)[18:23]


## ----TVCModel1_3, tidy=FALSE, echo=TRUE, message=FALSE, warning=FALSE----
M2 <- coxph(Surv(begin, end, event) ~ 
                qmv + qmvpostsea + qmvpostteu +
                coop + codec + eu9 + eu10 + eu12 + eu15 + thatcher + 
                agenda + backlog + 
                qmv_log + qmvpostsea_log + coop_log + 
                codec_log + thatcher_log + backlog_log, 
            data = GolubEUPData, ties = "efron") 


## ----TVCModel2_1, tidy=FALSE, message=FALSE------------------------------
Sim3 <- coxsimtvc(obj = M2, b = "qmv", btvc = "qmv_log",
                    qi = "First Difference", Xj = 1,
                    tfun = "log", from = 80, to = 2000, by = 10,
                    )


## ----TVCModel2_2, eval=FALSE, tidy=FALSE---------------------------------
## simGG(Sim3, xlab = "\nTime in Days", type = "ribbons", lsize = 0.5,
##       legend = FALSE, alpha = 0.3, nsim = 100)


## ----TVCModel1, tidy=FALSE, echo=FALSE, message=F, warning=F, fig.width=7, fig.height=4, out.width='0.95\\linewidth', dev='pdf'----
# Create first difference plots
simGG(Sim3, xlab = "\nTime in Days", alpha = 0.3,
      type = "ribbons", lsize = 0.5, legend = FALSE)


## ----TVCBacklog1, echo=TRUE, tidy=FALSE, eval=FALSE----------------------
## Sim4 <- coxsimtvc(obj = M2, b = "backlog", btvc = "backlog_log",
##                   qi = "Relative Hazard", Xj = seq(40, 200, 40),
##                   tfun = "log", from = 1200, to = 7000, by = 100,
##                   nsim = 100)
## 
## simGG(Sim4, xlab = "\nTime in Days", type = "ribbons",
##       leg.name = "Backlogged \n Items")


## ----TVCBacklog2, tidy=FALSE, echo=FALSE, message=F, warning=F, fig.width=7, fig.height=4, out.width='0.95\\linewidth', dev='pdf'----
# Create simtvc object for relative hazard
Sim4 <- coxsimtvc(obj = M2, b = "backlog", btvc = "backlog_log",
                  qi = "Relative Hazard", Xj = seq(40, 200, 40),
                  tfun = "log", from = 1200, to = 7000, by = 100,
                  nsim = 100) 

# Create relative hazard plot
simGG(Sim4, xlab = "\nTime in Days", type = "ribbons",
      leg.name = "Backlogged \n Items")


## ----FitKeeleModel, include=TRUE, tidy=FALSE,----------------------------
data("CarpenterFdaData")

M3 <- coxph(Surv(acttime, censor) ~  prevgenx + lethal + deathrt1 +
              acutediz + hosp01 + pspline(hospdisc, df = 4) + 
              pspline(hhosleng, df = 4) + mandiz01 + 
              femdiz01 + peddiz01 + orphdum + natreg + vandavg3 + 
              wpnoavg3 + pspline(condavg3, df = 4) + 
              pspline(orderent, df = 4) + pspline(stafcder, df = 4), 
            data = CarpenterFdaData)


## ----Spline1Fig1, eval=FALSE, tidy=FALSE---------------------------------
## XjFit <- seq(1100, 1700, by = 10)
## 
## XlFit <- setXl(Xj = XjFit, diff = 1)
## 
## Sim4 <- coxsimSpline(M3, bspline = "pspline(stafcder, df = 4)",
##                      bdata = CarpenterFdaData$stafcder,
##                      qi = "Hazard Ratio", nsim = 100,
##                      Xj = XjFit, Xl = XlFit)
## 
## simGG(Sim4, xlab = "\n Number of FDA Drug Review Staff",
##         title = "Central Interval\n", alpha = 0.1,
##         type = "lines")


## ----Spline1Fig2, echo=FALSE, tidy=FALSE, message=FALSE, warning=FALSE, fig.width=7, fig.height=4, out.width='0.95\\linewidth', cache=TRUE, dev='pdf'----
## Set Xj and Xl values
XjFit <- seq(1100, 1700, by = 10)
XlFit <- setXl(Xj = XjFit, diff = 1)

## Simulated Fitted Values
Sim5 <- coxsimSpline(M3, bspline = "pspline(stafcder, df = 4)", 
                     bdata = CarpenterFdaData$stafcder,
                     qi = "Hazard Ratio", nsim = 100,
                     Xj = XjFit, Xl = XlFit)

# Plot simulated values
Plot5_1 <- simGG(Sim5, xlab = "\n No. of FDA Drug Review Staff",
                  title = "Central Interval\n", alpha = 0.1, 
                  type = "lines")

Plot5_1 <- Plot5_1 + scale_y_continuous(breaks = c(0, 20, 40, 60), 
                                     limits = c(0, 60))

# Simulated Fitted Values: shortest probability interval
Sim6 <- coxsimSpline(M3, bspline = "pspline(stafcder, df = 4)", 
                     bdata = CarpenterFdaData$stafcder,
                     qi = "Hazard Ratio", nsim = 100,
                     Xj = XjFit, Xl = XlFit, 
                     spin = TRUE)

# Plot simulated values
Plot5_2 <- simGG(Sim6, xlab = "\n No. of FDA Drug Review Staff", ylab = "",
                title = "SPIn\n", alpha = 0.1, type = "lines")

# Place on the same scale as the central interval figure
Plot5_2 <- Plot5_2 + scale_y_continuous(breaks = c(0, 20, 40, 60), 
                                    limits = c(0, 60))

# Return combined plot
grid.arrange(Plot5_1, Plot5_2, ncol = 2)


## ----Spline2Fig1, eval=FALSE, tidy=FALSE---------------------------------
## library("ggplot2")
## 
## SimPlot + scale_y_continuous(breaks = c(0, 20, 40, 60), limits = c(0, 60))


