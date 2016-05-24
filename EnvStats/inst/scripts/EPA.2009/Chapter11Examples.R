#####
# File:     Script - EPA09, Chapter 11 Examples.R
#
# Purpose:  Reproduce Examples in Chapter 11 of the EPA Guidance Document
#
#           USEPA. (2009).  "Statistical Analysis of Ground Water 
#             Monitoring Data at RCRA Facilities, Unified Guidance."
#             EPA 530-R-09-007.
#             Office of Resource Conservation and Recovery, 
#             Program Information and Implementation Division.
#             March 2009.
#
# Author:   Steven P. Millard
#
# Last
# Updated:  September 3, 2012
#####

#NOTE:  Unless otherwise noted, differences between what is shown in the
#       EPA Guidance Document and results shown here are due to the 
#       EPA Guidance Document rounding intermediate results prior to 
#       the final result.

#######################################################################################


# Example 11-1, pp. 11-2 to 11-4
#-------------------------------

head(EPA.09.Ex.11.1.arsenic.df)
#  Arsenic.ppb Month Well
#1        22.9     1    1
#2         3.1     2    1
#3        35.7     3    1
#4         4.2     4    1
#5         2.0     1    2
#6         1.2     2    2

longToWide(EPA.09.Ex.11.1.arsenic.df, "Arsenic.ppb", "Month", "Well", 
	paste.row.name = TRUE, paste.col.name = TRUE)
#        Well.1 Well.2 Well.3 Well.4 Well.5 Well.6
#Month.1   22.9    2.0    2.0    7.8   24.9    0.3
#Month.2    3.1    1.2  109.4    9.3    1.3    4.8
#Month.3   35.7    7.8    4.5   25.9    0.8    2.8
#Month.4    4.2   52.0    2.5    2.0   27.0    1.2

summaryStats(Arsenic.ppb ~ Well, data = EPA.09.Ex.11.1.arsenic.df)
#  N   Mean      SD Median Min   Max
#1 4 16.475 15.7104  13.55 3.1  35.7
#2 4 15.750 24.3450   4.90 1.2  52.0
#3 4 29.600 53.2110   3.50 2.0 109.4
#4 4 11.250 10.2614   8.55 2.0  25.9
#5 4 13.500 14.4030  13.10 0.8  27.0
#6 4  2.275  1.9755   2.00 0.3   4.8

Resids <- lm(Arsenic.ppb ~ Well, data = EPA.09.Ex.11.1.arsenic.df)$residuals
new.df <- data.frame(EPA.09.Ex.11.1.arsenic.df, Resids)
Resids.mat <- longToWide(new.df, "Resids", "Month", "Well", 
	paste.row.name = TRUE, paste.col.name = TRUE)
Resids.mat 
#         Well.1 Well.2 Well.3 Well.4 Well.5 Well.6
#Month.1   6.425 -13.75  -27.6  -3.45   11.4 -1.975
#Month.2 -13.375 -14.55   79.8  -1.95  -12.2  2.525
#Month.3  19.225  -7.95  -25.1  14.65  -12.7  0.525
#Month.4 -12.275  36.25  -27.1  -9.25   13.5 -1.075

windows()
boxplot(Resids ~ Well, data = new.df, 
	names = paste("Well", levels(new.df$Well)),
	ylab = "Arsenic Residuals (ppb)",
	main = "Figure 11-1. Side-by-Side Box Plots of Arsenic Residuals")

rm(Resids, new.df, Resids.mat)


##################################################################################


# Example 11-2, pp. 11-7 to 11-8
#-------------------------------

varGroupTest(Arsenic.ppb ~ Well, data = EPA.09.Ex.11.1.arsenic.df)

#Results of Hypothesis Test
#--------------------------
#
#Null Hypothesis:                 Ratio of each pair of variances = 1
#
#Alternative Hypothesis:          At least one variance differs
#
#Test Name:                       Levene's Test for
#                                 Homogenity of Variance
#
#Estimated Parameter(s):          Well.1 =  246.8158
#                                 Well.2 =  592.6767
#                                 Well.3 = 2831.4067
#                                 Well.4 =  105.2967
#                                 Well.5 =  207.4467
#                                 Well.6 =    3.9025
#
#Data:                            Arsenic.ppb
#
#Grouping Variable:               Well
#
#Data Source:                     EPA.09.Ex.11.1.arsenic.df
#
#Sample Sizes:                    Well.1 = 4
#                                 Well.2 = 4
#                                 Well.3 = 4
#                                 Well.4 = 4
#                                 Well.5 = 4
#                                 Well.6 = 4
#
#Test Statistic:                  F = 4.564176
#
#Test Statistic Parameters:       num df   =  5
#                                 denom df = 18
#
#P-value:                         0.007294084


##################################################################################


# Example 11-3, pp. 11-9 to 11-10
#--------------------------------

Stats <- summaryStats(Arsenic.ppb ~ Well, data = EPA.09.Ex.11.1.arsenic.df, 
	combine.groups = F, digits = 3)
Stats
#  N   Mean     SD Median Min   Max
#1 4 16.480 15.710  13.55 3.1  35.7
#2 4 15.750 24.345   4.90 1.2  52.0
#3 4 29.600 53.211   3.50 2.0 109.4
#4 4 11.250 10.261   8.55 2.0  25.9
#5 4 13.500 14.403  13.10 0.8  27.0
#6 4  2.275  1.975   2.00 0.3   4.8


windows()
plot(Stats[, "Mean"], Stats[, "SD"], pch = 16,
	xlab = "Mean Arsenic (ppb)",
	ylab = "SD Arsenic (ppb)",
	main = "Figure 11-3. Arsenic Mean-Standard Deviation Plot")


Stats <- summaryStats(log(Arsenic.ppb) ~ Well, 
	data = EPA.09.Ex.11.1.arsenic.df, 
	combine.groups = F, digits = 3)
windows()
plot(Stats[, "Mean"], Stats[, "SD"], pch = 16,
	xlim = c(0, 3), ylim = c(0, 3),
	xlab = "Log-Mean Arsenic log(ppb)",
	ylab = "Log-SD Arsenic log(ppb)",
	main = "Figure 11-4. Log(Arsenic) Mean-Standard Deviation Plot")

rm(Stats)
