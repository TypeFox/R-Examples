#############
# Time-interactive effect examples from JSS paper
# Christopher Gandrud
# Updated 6 April 2014
#############

# Load packages
library("survival")
library("simPH")
library("ggplot2")
library("gridExtra")


##### Illustration of time-varying interactive effects ######
# Load Golub & Steunenberg (2007) data. The data is included with simPH.
data("GolubEUPData")

# Examine the data's format
head(GolubEUPData[, 2:5])

# Expand data into equally spaced time intervals
 GolubEUPData <- SurvExpand(GolubEUPData, GroupVar = 'caseno',
                      Time = 'begin', Time2 = 'end', event = 'event') 

# Examine the data's new format
head(GolubEUPData[, 1:4])

# Create natural log-time interactions
BaseVars <- c('qmv', 'backlog', 'coop', 'codec', 'qmvpostsea', 'thatcher')
GolubEUPData <- tvc(GolubEUPData, b = BaseVars, tvar = 'end', tfun = 'log')

# Estimate model
M2 <- coxph(Surv(begin, end, event) ~ qmv + qmvpostsea + qmvpostteu +
                coop + codec + eu9 + eu10 + eu12 + eu15 + thatcher +
                agenda + backlog + qmv_log + qmvpostsea_log + coop_log +
                codec_log + thatcher_log + backlog_log,
            data = GolubEUPData, ties = "efron")

## Create simtvc object for first difference (central interval)
Sim3 <- coxsimtvc(obj = M2, b = "qmv", btvc = "qmv_log",
                    qi = "First Difference", Xj = 1,
                    tfun = "log", from = 80, to = 2000,
                    by = 10, ci = 0.95)

# Create first difference plots
simGG(Sim3, xlab = "\nTime in Days", 
       title = "Central Interval\n", alpha = 0.3,
       type = "ribbons", lsize = 0.5, legend = FALSE)

# Create simtvc object for relative hazard
Sim4 <- coxsimtvc(obj = M2, b = "backlog", btvc = "backlog_log",
                  qi = "Relative Hazard", Xj = seq(40, 200, 40),
                  tfun = "log", from = 1200, to = 7000, by = 100,
                  nsim = 200)

# Create relative hazard plot
simGG(Sim4, xlab = "\nTime in Days", type = "ribbons",
      leg.name = "Backlogged \n Items")
