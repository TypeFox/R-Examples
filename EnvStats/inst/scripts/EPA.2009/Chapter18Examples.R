#####
# File:     Script - EPA09, Chapter 18 Examples.R
#
# Purpose:  Reproduce Examples in Chapter 18 of the EPA Guidance Document
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

# Example 18-1, pp. 18-9 to 18-10
#--------------------------------

EPA.09.Ex.18.1.arsenic.df
#   Year Sampling.Period Arsenic.ppb
#1     1      Background        12.6
#2     1      Background        30.8
#3     1      Background        52.0
#4     1      Background        28.1
#5     2      Background        33.3
#6     2      Background        44.0
#7     2      Background         3.0
#8     2      Background        12.8
#9     3      Background        58.1
#10    3      Background        12.6
#11    3      Background        17.6
#12    3      Background        25.3
#13    4      Compliance        48.0
#14    4      Compliance        30.3
#15    4      Compliance        42.5
#16    4      Compliance        15.0


# Summary statistics
#-------------------

summaryStats(Arsenic.ppb ~ Sampling.Period, 
	data = EPA.09.Ex.18.1.arsenic.df,
	digits = 2, combine.groups = TRUE, 
	stats.in.rows = TRUE)
#
#       Background Compliance Combined
#N      12          4         16      
#Mean   27.52      33.95      29.12   
#SD     17.1       14.64      16.3    
#Median 26.7       36.4       29.2    
#Min     3         15          3      
#Max    58.1       48         58.1


# Step 1:  Shapiro-Wilk Test on Background Data
#----------------------------------------------

gofTest(Arsenic.ppb ~ 1, data = EPA.09.Ex.18.1.arsenic.df,
	subset = Sampling.Period == "Background")
#
#Results of Goodness-of-Fit Test
#-------------------------------
#
#Test Method:                     Shapiro-Wilk GOF
#
#Hypothesized Distribution:       Normal
#
#Estimated Parameter(s):          mean = 27.51667
#                                 sd   = 17.10119
#
#Estimation Method:               mvue
#
#Data:                            Arsenic.ppb
#
#Subset With:                     Sampling.Period == "Background"
#
#Data Source:                     EPA.09.Ex.18.1.arsenic.df
#
#Sample Size:                     12
#
#Test Statistic:                  W = 0.94695
#
#Test Statistic Parameter:        n = 12
#
#P-value:                         0.5929102
#
#Alternative Hypothesis:          True cdf does not equal the
#                                 Normal Distribution.


# Step 2 - Predicition Limit for 4 future observations
#-----------------------------------------------------

Sampling.Period <- EPA.09.Ex.18.1.arsenic.df$Sampling.Period
As              <- EPA.09.Ex.18.1.arsenic.df$Arsenic.ppb

pred.int.list <- 
	predIntNorm(As[Sampling.Period == "Background"], 
	k = 4, pi.type = "upper", conf.level = 0.95) 

pred.int.list
#
#Results of Distribution Parameter Estimation
#--------------------------------------------
#
#Assumed Distribution:            Normal
#
#Estimated Parameter(s):          mean = 27.51667
#                                 sd   = 17.10119
#
#Estimation Method:               mvue
#
#Data:                            As[Sampling.Period == "Background"]
#
#Sample Size:                     12
#
#Prediction Interval Method:      Bonferroni
#
#Prediction Interval Type:        upper
#
#Confidence Level:                95%
#
#Number of Future Observations:   4
#
#Prediction Interval:             LPL =     -Inf
#                                 UPL = 73.67237

names(pred.int.list)
#[1] "distribution" "sample.size"  "parameters"  
#[4] "n.param.est"  "method"       "data.name"   
#[7] "bad.obs"      "interval"

pred.int.list$interval$limits["UPL"]
#     UPL 
#73.67237

any(As[Sampling.Period == "Compliance"] > 
	pred.int.list$interval$limits["UPL"])
#[1] FALSE

   
rm(Sampling.Period, As, pred.int.list)


#######################################################################################

# Figure 18-1, p. 18-12
#----------------------

windows()
par(mfrow = c(2,2), mar = c(3, 3, 2, 1), 
	mgp = c(1.5, 0.5, 0), oma = c(0, 0, 3, 0))

# 2 future values (p = 2)
#------------------------
plotPredIntNormTestPowerCurve(n=8, k=1, n.mean = 2,
	range.delta.over.sigma=c(0, 6), 
	pi.type="upper", conf.level=0.99,
	plot.lty = 1, plot.col = "black",
	ylim = c(0, 1), 
	xlab = "SD Units Over BG", ylab = "Power", 
	main = "p = 2")

plotPredIntNormSimultaneousTestPowerCurve(n=8, k=1, m=2, 
	range.delta.over.sigma=c(0, 6), 
	pi.type="upper", conf.level=0.99, 
	add=T, plot.lty = 2, plot.col = "red")

plotPredIntNormTestPowerCurve(n=8, k=2, n.mean = 1,
	range.delta.over.sigma=c(0, 6), 
	pi.type="upper", conf.level=0.99,
	add = T, plot.lty = 3, plot.col = "blue")

legend(2.25, 0.35, 
	c("PL on mean", "PL w/ retests (1-of-2)", "PL on values"), 
	lty = 1:3, col = c("black", "red", "blue"),
	cex = 0.85, lwd = 1.5, bty = "n")


# 3 future values (p = 3)
#------------------------

plotPredIntNormTestPowerCurve(n=8, k=1, n.mean = 3,
	range.delta.over.sigma=c(0, 6), 
	pi.type="upper", conf.level=0.99,
	plot.lty = 1, plot.col = "black",
	ylim = c(0, 1), 
	xlab = "SD Units Over BG", ylab = "Power", 
	main = "p = 3")

plotPredIntNormSimultaneousTestPowerCurve(n=8, k=1, m=3, 
	range.delta.over.sigma=c(0, 6), 
	pi.type="upper", conf.level=0.99, 
	add=T, plot.lty=2, plot.col = "red")

plotPredIntNormTestPowerCurve(n=8, k=3, n.mean = 1,
	range.delta.over.sigma=c(0, 6), 
	pi.type="upper", conf.level=0.99,
	add = T, plot.lty = 3, plot.col = "blue")

legend(2.25, 0.35, 
	c("PL on mean", "PL w/ retests (1-of-3)", "PL on values"), 
	lty = 1:3, col = c("black", "red", "blue"),
	cex = 0.85, lwd = 1.5, bty = "n")


# 4 future values (p = 4)
#------------------------

plotPredIntNormTestPowerCurve(n=8, k=1, n.mean = 4,
	range.delta.over.sigma=c(0, 6), 
	pi.type="upper", conf.level=0.99,
	plot.lty = 1, plot.col = "black",
	ylim = c(0, 1), 
	xlab = "SD Units Over BG", ylab = "Power", 
	main = "p = 4")

plotPredIntNormSimultaneousTestPowerCurve(n=8, k=1, m=4, 
	range.delta.over.sigma=c(0, 6), 
	pi.type="upper", conf.level=0.99, 
	add=T, plot.lty = 2, plot.col = "red")

plotPredIntNormTestPowerCurve(n=8, k=4, n.mean = 1,
	range.delta.over.sigma=c(0, 6), 
	pi.type="upper", conf.level=0.99,
	add = T, plot.lty = 3, plot.col = "blue")

legend(2.25, 0.35, 
	c("PL on mean", "PL w/ retests (1-of-4)", "PL on values"), 
	lty = 1:3, col = c("black", "red", "blue"),
	cex = 0.85, lwd = 1.5, bty = "n")

mtext("Figure 18-1. Comparison of Prediction Limits", 
	outer = TRUE, cex = 1.25, line = 2)
mtext(expression(paste("(BG = 8, ", alpha, " = 0.01, 1 test)")),
	outer = TRUE, cex = 1.25, line = 0.5)



# Figure 18-2, p. 18-13
#----------------------

windows()
par(mfrow = c(2,2), mar = c(3, 3, 2, 1), 
	mgp = c(1.5, 0.5, 0), oma = c(0, 0, 3, 0))

# 2 future values (p = 2)
#------------------------
plotPredIntNormTestPowerCurve(n=20, k=1, n.mean = 2,
	range.delta.over.sigma=c(0, 5), 
	pi.type="upper", conf.level=0.95,
	plot.lty = 1, plot.col = "black",
	ylim = c(0, 1), 
	xlab = "SD Units Over BG", ylab = "Power", 
	main = "p = 2")

plotPredIntNormSimultaneousTestPowerCurve(n=20, k=1, m=2, 
	range.delta.over.sigma=c(0, 5), 
	pi.type="upper", conf.level=0.95, 
	add=T, plot.col = "red", plot.lty=3)

plotPredIntNormTestPowerCurve(n=20, k=2, n.mean = 1,
	range.delta.over.sigma=c(0, 5), 
	pi.type="upper", conf.level=0.95,
	add = T, plot.lty = 2, plot.col = "blue")

legend(1.75, 0.35, 
	c("PL on mean", "PL w/ retests (1-of-2)", "PL on values"), 
	lty = c(1, 3, 2), col = c("black", "red", "blue"),
	cex = 0.85, lwd = 1.5, bty = "n")


# 3 future values (p = 3)
#------------------------

plotPredIntNormTestPowerCurve(n=20, k=1, n.mean = 3,
	range.delta.over.sigma=c(0, 5), 
	pi.type="upper", conf.level=0.95,
	plot.lty = 1, plot.col = "black",
	ylim = c(0, 1), 
	xlab = "SD Units Over BG", ylab = "Power", 
	main = "p = 3")

plotPredIntNormSimultaneousTestPowerCurve(n=20, k=1, m=3, 
	range.delta.over.sigma=c(0, 5), 
	pi.type="upper", conf.level=0.95, 
	add=T, plot.col = "red", plot.lty=3)

plotPredIntNormTestPowerCurve(n=20, k=3, n.mean = 1,
	range.delta.over.sigma=c(0, 5), 
	pi.type="upper", conf.level=0.95,
	add = T, plot.lty = 2, plot.col = "blue")

legend(1.75, 0.35, 
	c("PL on mean", "PL w/ retests (1-of-3)", "PL on values"), 
	lty = c(1, 3, 2), col = c("black", "red", "blue"),
	cex = 0.85, lwd = 1.5, bty = "n")


# 4 future values (p = 4)
#------------------------

plotPredIntNormTestPowerCurve(n=20, k=1, n.mean = 4,
	range.delta.over.sigma=c(0, 5), 
	pi.type="upper", conf.level=0.95,
	plot.lty = 1, plot.col = "black",
	ylim = c(0, 1), 
	xlab = "SD Units Over BG", ylab = "Power", 
	main = "p = 4")

plotPredIntNormSimultaneousTestPowerCurve(n=20, k=1, m=4, 
	range.delta.over.sigma=c(0, 5), 
	pi.type="upper", conf.level=0.95, 
	add=T, plot.col = "red", plot.lty=3)

plotPredIntNormTestPowerCurve(n=20, k=4, n.mean = 1,
	range.delta.over.sigma=c(0, 5), 
	pi.type="upper", conf.level=0.95,
	add = T, plot.lty = 2, plot.col = "blue")

legend(1.75, 0.35, 
	c("PL on mean", "PL w/ retests (1-of-4)", "PL on values"), 
	lty = c(1, 3, 2), col = c("black", "red", "blue"),
	cex = 0.85, lwd = 1.5, bty = "n")

mtext("Figure 18-2. Comparison of Prediction Limits", 
	outer = TRUE, cex = 1.25, line = 2)
mtext(expression(paste("(BG = 20, ", alpha, " = 0.05, 1 test)")),
	outer = TRUE, cex = 1.25, line = 0.5)


#####################################################################################################

# Example 18-2, pp. 18-15 to 18-16
#---------------------------------

# Raw data and summary statistics
#--------------------------------

EPA.09.Ex.18.2.chrysene.df
#   Month   Well  Well.type Chrysene.ppb
#1      1 Well.1 Background          6.9
#2      2 Well.1 Background         27.3
#3      3 Well.1 Background         10.8
#4      4 Well.1 Background          8.9
#5      1 Well.2 Background         15.1
#6      2 Well.2 Background          7.2
#7      3 Well.2 Background         48.4
#8      4 Well.2 Background          7.8
#9      1 Well.3 Compliance         68.0
#10     2 Well.3 Compliance         48.9
#11     3 Well.3 Compliance         30.1
#12     4 Well.3 Compliance         38.1

longToWide(EPA.09.Ex.18.2.chrysene.df, 
	"Chrysene.ppb", "Month", "Well",
	paste.row.name = TRUE)
#
#        Well.1 Well.2 Well.3
#Month.1    6.9   15.1   68.0
#Month.2   27.3    7.2   48.9
#Month.3   10.8   48.4   30.1
#Month.4    8.9    7.8   38.1


summaryStats(Chrysene.ppb ~ Well, 
	data = EPA.09.Ex.18.2.chrysene.df,
	combine.groups = FALSE, 
	digits = 2, stats.in.rows = TRUE)
#
#       Well.1 Well.2 Well.3
#N       4      4      4    
#Mean   13.48  19.62  46.27 
#SD      9.35  19.52  16.4  
#Median  9.85  11.45  43.5  
#Min     6.9    7.2   30.1  
#Max    27.3   48.4   68     

summaryStats(Chrysene.ppb ~ Well.type,
	data = EPA.09.Ex.18.2.chrysene.df, 
	combine.groups = FALSE, digits = 2, 
	stats.in.rows = TRUE)
#
#       Background Compliance
#N       8          4        
#Mean   16.55      46.27     
#SD     14.54      16.4      
#Median  9.85      43.5      
#Min     6.9       30.1      
#Max    48.4       68      
#Max    48.4       68


# Step 1 - Check Distribution Assumptions for 
#          Background Observations.
#--------------------------------------------


# Shapiro-Wilk Test for Normal Distribution
#-------------------------------------------

gofTest(Chrysene.ppb ~ 1, data = EPA.09.Ex.18.2.chrysene.df,
	subset = Well.type == "Background")
#
#
#Results of Goodness-of-Fit Test
#-------------------------------
#
#Test Method:                     Shapiro-Wilk GOF
#
#Hypothesized Distribution:       Normal
#
#Estimated Parameter(s):          mean = 16.55000
#                                 sd   = 14.54441
#
#Estimation Method:               mvue
#
#Data:                            Chrysene.ppb
#
#Subset With:                     Well.type == "Background"
#
#Data Source:                     EPA.09.Ex.18.2.chrysene.df
#
#Sample Size:                     8
#
#Test Statistic:                  W = 0.7289006
#
#Test Statistic Parameter:        n = 8
#
#P-value:                         0.004759859
#
#Alternative Hypothesis:          True cdf does not equal the
#                                 Normal Distribution.


# Shapiro-Wilk Test for Lognormal distribution
#---------------------------------------------

gofTest(Chrysene.ppb ~ 1, data = EPA.09.Ex.18.2.chrysene.df,
	subset = Well.type == "Background", dist = "lnorm")

#
#Results of Goodness-of-Fit Test
#-------------------------------
#
#Test Method:                     Shapiro-Wilk GOF
#
#Hypothesized Distribution:       Lognormal
#
#Estimated Parameter(s):          meanlog = 2.5533006
#                                 sdlog   = 0.7060038
#
#Estimation Method:               mvue
#
#Data:                            Chrysene.ppb
#
#Subset With:                     Well.type == "Background"
#
#Data Source:                     EPA.09.Ex.18.2.chrysene.df
#
#Sample Size:                     8
#
#Test Statistic:                  W = 0.8546352
#
#Test Statistic Parameter:        n = 8
#
#P-value:                         0.1061057
#
#Alternative Hypothesis:          True cdf does not equal the
#                                 Lognormal Distribution.


# Raw data and summary statistics for log-transformed data
#---------------------------------------------------------

summaryStats(log(Chrysene.ppb) ~ Well, 
	data = EPA.09.Ex.18.2.chrysene.df,
	combine.groups = FALSE, 
	digits = 3, stats.in.rows = TRUE)
#
#       Well.1 Well.2 Well.3
#N      4      4      4     
#Mean   2.451  2.656  3.789 
#SD     0.599  0.881  0.349 
#Median 2.283  2.384  3.765 
#Min    1.932  1.974  3.405 
#Max    3.307  3.879  4.22  
 

summaryStats(log(Chrysene.ppb) ~ Well.type, 
	data = EPA.09.Ex.18.2.chrysene.df,
	combine.groups = FALSE, 
	digits = 3, stats.in.rows = TRUE)
#
#       Background Compliance
#N      8          4         
#Mean   2.553      3.789     
#SD     0.706      0.349     
#Median 2.283      3.765     
#Min    1.932      3.405     
#Max    3.879      4.22


# Compute upper 99% Prediction Interval for mean of 4 future observations
# based on log-transformed Chrysene data for Background Wells
#------------------------------------------------------------------------

Well.type <- EPA.09.Ex.18.2.chrysene.df$Well.type
Chrysene  <- EPA.09.Ex.18.2.chrysene.df$Chrysene.ppb

predIntNorm(log(Chrysene)[Well.type == "Background"], k = 1, n.mean = 4, 
	method = "exact", pi.type = "upper", conf.level = 0.99)
#
#Results of Distribution Parameter Estimation
#--------------------------------------------
#
#Assumed Distribution:            Normal
#
#Estimated Parameter(s):          mean = 2.5533006
#                                 sd   = 0.7060038
#
#Estimation Method:               mvue
#
#Data:                            log(Chrysene)[Well.type == "Background"]
#
#Sample Size:                     8
#
#Prediction Interval Method:      exact
#
#Prediction Interval Type:        upper
#
#Confidence Level:                99%
#
#Number of Future Averages:       1
#
#Sample Size for Averages:        4
#
#Prediction Interval:             LPL =     -Inf
#                                 UPL = 3.849427

# Compare 3.85 to the observed mean of the log-transformed Compliance Well
# observations:  3.789 is less than 3.85, so no evidence of contamination.



#NOTE:  You can also compute the prediction limit for the *geometric* mean 
#       of the next 4 observations for the untransformed data:      

predIntLnorm(Chrysene[Well.type == "Background"], k = 1, n.geomean = 4, 
	method = "exact", pi.type = "upper", conf.level = 0.99)
#
#Results of Distribution Parameter Estimation
#--------------------------------------------
#
#Assumed Distribution:            Lognormal
#
#Estimated Parameter(s):          meanlog = 2.5533006
#                                 sdlog   = 0.7060038
#
#Estimation Method:               mvue
#
#Data:                            Chrysene[Well.type == "Background"]
#
#Sample Size:                     8
#
#Prediction Interval Method:      exact
#
#Prediction Interval Type:        upper
#
#Confidence Level:                99%
#
#Number of Future Averages:       1
#
#Sample Size for Averages:        4
#
#Prediction Interval:             LPL =  0.00000
#                                 UPL = 46.96613

# Compare 46.97 to the observed geometric mean of the Compliance Well
# observations:  

geo.mean(Chrysene[Well.type == "Compliance"])
#[1] 44.19034

#44.19 is less than 46.97, so no evidence of contamination.


rm(Well.type, Chrysene)


#####################################################################################################

# Example 18-3, p. 18-19
#-----------------------

head(EPA.09.Ex.18.3.TCE.df)
#  Month Well  Well.type TCE.ppb.orig TCE.ppb Censored
#1     1 BW-1 Background           <5       5     TRUE
#2     2 BW-1 Background           <5       5     TRUE
#3     3 BW-1 Background            8       8    FALSE
#4     4 BW-1 Background           <5       5     TRUE
#5     5 BW-1 Background            9       9    FALSE
#6     6 BW-1 Background           10      10    FALSE


longToWide(EPA.09.Ex.18.3.TCE.df, "TCE.ppb.orig",
	"Month", "Well", paste.row.name = TRUE)
#
#        BW-1 BW-2 BW-3 CW-4
#Month.1   <5    7   <5     
#Month.2   <5  6.5   <5     
#Month.3    8   <5 10.5  7.5
#Month.4   <5    6   <5   <5
#Month.5    9   12   <5    8
#Month.6   10   <5    9   14



# Construct a nonparametric prediction interval for 
# the next 4 future observations using the maximum value of the
# background observations as the upper prediction limit.

Well.type <- EPA.09.Ex.18.3.TCE.df$Well.type
TCE       <- EPA.09.Ex.18.3.TCE.df$TCE.ppb

predIntNpar(TCE[Well.type == "Background"], m = 4, 
	pi.type = "upper", lb = 0) 
#
#Results of Distribution Parameter Estimation
#--------------------------------------------
#
#Assumed Distribution:            None
#
#Data:                            TCE[Well.type == "Background"]
#
#Sample Size:                     18
#
#Prediction Interval Method:      Exact
#
#Prediction Interval Type:        upper
#
#Confidence Level:                82%
#
#Prediction Limit Rank(s):        18 
#
#Number of Future Observations:   4
#
#Prediction Interval:             LPL =  0
#                                 UPL = 12

# Month 6 TCE measure at compliance well is 14 ppb > 12, so conclude
# there is evidence of contamination at the 100%-82% = 18% Type I Error level.


rm(Well.type, TCE)

#######################################################################################################

# Example 18-4, p. 18-21 to 18-22
#--------------------------------

head(EPA.09.Ex.18.4.xylene.df)
#  Month   Well  Well.type Xylene.ppb.orig Xylene.ppb Censored
#1     1 Well.1 Background              <5        5.0     TRUE
#2     2 Well.1 Background              <5        5.0     TRUE
#3     3 Well.1 Background             7.5        7.5    FALSE
#4     4 Well.1 Background              <5        5.0     TRUE
#5     5 Well.1 Background              <5        5.0     TRUE
#6     6 Well.1 Background              <5        5.0     TRUE
   

longToWide(EPA.09.Ex.18.4.xylene.df, "Xylene.ppb.orig",
	"Month", "Well", paste.row.name = TRUE)
#        Well.1 Well.2 Well.3 Well.4
#Month.1     <5    9.2     <5       
#Month.2     <5     <5    5.4       
#Month.3    7.5     <5    6.7       
#Month.4     <5    6.1     <5       
#Month.5     <5      8     <5       
#Month.6     <5    5.9     <5     <5
#Month.7    6.4     <5     <5    7.8
#Month.8      6     <5     <5   10.4

 

# Construct a nonparametric prediction interval for 
# the median of the next 3 future observations 
# using the maximum value of the background observations 
# as the upper prediction limit.
#
# This is equivalent to constructing a nonparametric 
# prediction interval that must hold at least 2 of the next 3
# future observations.


Well.type   <- EPA.09.Ex.18.4.xylene.df$Well.type
Xylene      <- EPA.09.Ex.18.4.xylene.df$Xylene.ppb

predIntNparSimultaneous(Xylene[Well.type == "Background"], 
	k = 2, m = 3, pi.type = "upper", lb = 0) 
#
#Results of Distribution Parameter Estimation
#--------------------------------------------
#
#Assumed Distribution:            None
#
#Data:                            Xylene[Well.type == "Background"]
#
#Sample Size:                     24
#
#Prediction Interval Method:      exact 
#
#Prediction Interval Type:        upper
#
#Confidence Level:                 99%
#
#Prediction Limit Rank(s):        24 
#
#Minimum Number of
#Future Observations
#Interval Should Contain:         2
#
#Total Number of
#Future Observations:             3
#
#Prediction Interval:             LPL = 0.0
#                                 UPL = 9.2

EPA.09.Ex.18.4.xylene.df[Well.type == "Compliance", 
	c("Month", "Well", "Xylene.ppb.orig")]
#   Month   Well Xylene.ppb.orig
#25     1 Well.4                
#26     2 Well.4                
#27     3 Well.4                
#28     4 Well.4                
#29     5 Well.4                
#30     6 Well.4              <5
#31     7 Well.4             7.8
#32     8 Well.4            10.4 

# The Month 8 observation at the Complance well is 10.4 ppb of Xylene, 
# which is greater than the upper prediction limit of 9.2 ppb, so
# conclude there is evidence of contamination at the 
# 100% - 99% = 1% Type I Error Level


rm(Well.type, Xylene)

############################################################################

