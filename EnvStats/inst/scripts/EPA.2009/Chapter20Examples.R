#####
# File:     Script - EPA09, Chapter 20 Examples.R
#
# Purpose:  Reproduce Examples in Chapter 20 of the EPA Guidance Document
#
#           USEPA. (2009).  "Statistical Analysis of Ground Water 
#             Monitoring Data at RCRA Facilities, Unified Guidance."
#             EPA 530-R-09-007.
#             Office of Resource Conservation and Recovery, 
#             Program Information and Implementation Division.
#             March 2009.
#
#           Errata are listed in:
#           USEPA. (2010).  "Errata Sheet - March 2009 Unified Guidance."
#             EPA 530/R-09-007a.
#             Office of Resource Conservation and Recovery, 
#             Program Information and Implementation Division.
#             August 9, 2010.
#
# Author:   Steven P. Millard
#
# Last
# Updated:  October 3, 2012
#####

#NOTE:  Unless otherwise noted, differences between what is shown in the
#       EPA Guidance Document and results shown here are due to the 
#       EPA Guidance Document rounding intermediate results prior to 
#       the final result.

#######################################################################################

library(EnvStats)

#NOTE:  Functions for producing quality control statistics and graphs
#       are available in the R package qcc.  You must have this 
#       package installed in order to load it.

library(qcc)

#######################################################################################

# Example 20-1, pp. 20-4 to 20-6
#-------------------------------

head(EPA.09.Ex.20.1.nickel.df)
#  Month Year   Period Nickel.ppb
#1     1 1995 Baseline       32.8
#2     2 1995 Baseline       15.2
#3     3 1995 Baseline       13.5
#4     4 1995 Baseline       39.6
#5     5 1995 Baseline       37.1
#6     6 1995 Baseline       10.4

longToWide(EPA.09.Ex.20.1.nickel.df, "Nickel.ppb", 
	"Month", "Period", paste.row.name = TRUE)
#        Baseline Compliance
#Month.1     32.8       19.0
#Month.2     15.2       34.5
#Month.3     13.5       17.8
#Month.4     39.6       23.6
#Month.5     37.1       34.8
#Month.6     10.4       28.8
#Month.7     31.9       23.7
#Month.8     20.6       81.8


# Summary statistics by Period
#-----------------------------

summaryStats(Nickel.ppb ~ Period, 
	data = EPA.09.Ex.20.1.nickel.df,
	combine.groups = TRUE, stats.in.rows = TRUE, 
	digits = 1)
#       Baseline Compliance Combined
#N       8        8         16      
#Mean   25.1     33         29.1    
#SD     11.5     20.7       16.7    
#Median 26.2     26.2       26.2    
#Min    10.4     17.8       10.4    
#Max    39.6     81.8       81.8  


# Step 1:  Shapiro-Wilk Test on Background Data
#----------------------------------------------

gof.list <- gofTest(Nickel.ppb ~ 1, 
	data = EPA.09.Ex.20.1.nickel.df, 
	subset = Period == "Baseline")
gof.list
#
#Results of Goodness-of-Fit Test
#-------------------------------
#
#Test Method:                     Shapiro-Wilk GOF
#
#Hypothesized Distribution:       Normal
#
#Estimated Parameter(s):          mean = 25.13750
#                                 sd   = 11.51818
#
#Estimation Method:               mvue
#
#Data:                            Nickel.ppb
#
#Subset With:                     Period == "Baseline"
#
#Data Source:                     EPA.09.Ex.20.1.nickel.df
#
#Sample Size:                     8
#
#Test Statistic:                  W = 0.8964134
#
#Test Statistic Parameter:        n = 8
#
#P-value:                         0.2681304
#
#Alternative Hypothesis:          True cdf does not equal the
#                                 Normal Distribution.


windows()
plot(gof.list, plot.type = "Q-Q Plot",
	ylab = "Nickel Concentration (ppb)",
	main = "Figure 20-1. Probability Plot of Baseline Nickel Data")



# Steps 2-6 Create Shewhart and CUSUM control charts 
#---------------------------------------------------

Month <- EPA.09.Ex.20.1.nickel.df$Month
Period <- EPA.09.Ex.20.1.nickel.df$Period
Nickel <- EPA.09.Ex.20.1.nickel.df$Nickel.ppb

Nickel.Baseline   <- Nickel[Period == "Baseline"]
Nickel.Compliance <- Nickel[Period == "Compliance"]
Baseline.mean <- mean(Nickel.Baseline)
Baseline.sd   <- sd(Nickel.Baseline)
h <- 5


# CUSUM chart
#------------
windows()
cusum.Nickel <- cusum(data = Nickel.Compliance, 
	sizes = 1, center = Baseline.mean, 
	std.dev = Baseline.sd, 
	decision.interval = h, se.shift = 1,
	data.name = "Nickel Conc (ppb)", 
	xlab = "Month",
	title = "Figure 20-2a. Cusum Control Chart for Nickel Measurements") 


# Combined Shewhart-CUSUM chart
#------------------------------

cusum.obs   <- Baseline.mean + Baseline.sd * cusum.Nickel$pos

windows()
par(mar = c(5, 4, 4, 4) + 0.1)
plot(1:length(Nickel.Compliance), Nickel.Compliance, 
	type = "n", lty = 2, ylim = c(0, 100),
	xlab = "Sampling Event", ylab = "Nickel Conc (ppb)",
	main = "Figure 20-2. Shewhart-CUSUM Control Chart for Nickel Measurements")
points(1:length(Nickel.Compliance), Nickel.Compliance, 
	type = "b", lty = 2, lwd = 2)
points(1:length(Nickel.Compliance), cusum.obs, 
	type = "b", pch = 15, col = 2, lwd = 2)
abline(h = Baseline.mean + h * Baseline.sd, col = 4, lwd = 3)
axis(side = 4, at = Baseline.mean + h * Baseline.sd, labels = "h-limit", 
	tick = FALSE, las = 2)
legend(5, 15, c("Nickel", "CUSUM"), pch = c(1, 15), 
	lty = 2:1, col = 1:2, lwd = 2, bty = "n")

rm(Month, Period, Nickel, Nickel.Baseline, gof.list, 
	Nickel.Compliance, Baseline.mean, Baseline.sd, 
	h, cusum.Nickel, cusum.obs)
