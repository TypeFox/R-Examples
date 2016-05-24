#############
# Penalized spline exampes from JSS paper
# Christopher Gandrud
# Updated 5 April 2014
#############

# Load packages
library("survival")
library("simPH")
library("ggplot2")
library("gridExtra")

# Load Carpenter (2002) data. The data is included with simPH.
data("CarpenterFdaData")

# Run basic model
# From Keele (2010) replication source code. Used to create Table 7.
M3 <- coxph(Surv(acttime, censor) ~  prevgenx + lethal + deathrt1 +
              acutediz + hosp01 + pspline(hospdisc, df = 4) + 
              pspline(hhosleng, df = 4) + mandiz01 + 
              femdiz01 + peddiz01 + orphdum + natreg + vandavg3 + 
              wpnoavg3 + pspline(condavg3, df = 4) + 
              pspline(orderent, df = 4) + pspline(stafcder, df = 4), 
              data = CarpenterFdaData)

## Set Xj and Xl values
XjFit <- seq(1100, 1700, by = 10)
XlFit <- setXl(Xj = XjFit, diff = 1)

## Simulated Fitted Values
Sim5 <- coxsimSpline(M3, bspline = "pspline(stafcder, df = 4)", 
                     bdata = CarpenterFdaData$stafcder,
                     qi = "Hazard Ratio",
                     Xj = XjFit, Xl = XlFit)

# Plot simulated values
Plot5_1 <- simGG(Sim5, xlab = "\n Number of FDA Drug Review Staff",
                  title = "Central Interval\n", alpha = 0.1, 
                  type = "lines")

Plot5_1 <- Plot5_1 + scale_y_continuous(breaks = c(0, 20, 40, 60), 
                                     limits = c(0, 60))

# Simulated Fitted Values: shortest probability interval
Sim6 <- coxsimSpline(M3, bspline = "pspline(stafcder, df = 4)", 
                     bdata = CarpenterFdaData$stafcder,
                     qi = "Hazard Ratio",
                     Xj = XjFit, Xl = XlFit, 
                     spin = TRUE)

# Plot simulated values
Plot5_2 <- simGG(Sim6, xlab = "\n Number of FDA Drug Review Staff", ylab = "",
                title = "SPIn\n", alpha = 0.1, type = "lines")

# Place on the same scale as the central interval figure
Plot5_2 <- Plot5_2 + scale_y_continuous(breaks = c(0, 20, 40, 60), 
                                    limits = c(0, 60))

# Return combined plot
grid.arrange(Plot5_1, Plot5_2, ncol = 2)