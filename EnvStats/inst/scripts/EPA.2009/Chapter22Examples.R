#####
# File:     Script - EPA09, Chapter 22 Examples.R
#
# Purpose:  Reproduce Examples in Chapter 22 of the EPA Guidance Document
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
# Updated:  October 4, 2012
#####

#NOTE:  Unless otherwise noted, differences between what is shown in the
#       EPA Guidance Document and results shown here are due to the 
#       EPA Guidance Document rounding intermediate results prior to 
#       the final result.

#######################################################################################

library(EnvStats)

# Example 22-1, pp. 22-6 to 22-8
#-------------------------------

EPA.09.Ex.22.1.VC.df
#   Year Quarter     Period Well VC.ppb
#1     1       1 Background GW-1    6.3
#2     1       2 Background GW-1    9.5
#3     1       3 Background GW-1    8.1
#4     1       4 Background GW-1   11.9
#5     2       1 Compliance GW-1    7.3
#6     2       2 Compliance GW-1   11.2
#7     2       3 Compliance GW-1    6.0
#8     2       4 Compliance GW-1    7.5
#9     3       1 Compliance GW-1    8.4
#10    3       2 Compliance GW-1    6.4
#11    3       3 Compliance GW-1    8.9
#12    3       4 Compliance GW-1    4.9
#13    4       1 Compliance GW-1    9.6
#14    4       2 Compliance GW-1    9.7
#15    4       3 Compliance GW-1    8.7
#16    4       4 Compliance GW-1    8.7
#17    1       1 Background GW-2    5.9
#18    1       2 Background GW-2    3.0
#19    1       3 Background GW-2    8.8
#20    1       4 Background GW-2   12.0
#21    2       1 Compliance GW-2   11.2
#22    2       2 Compliance GW-2    8.6
#23    2       3 Compliance GW-2   12.6
#24    2       4 Compliance GW-2    7.2
#25    3       1 Compliance GW-2   13.8
#26    3       2 Compliance GW-2    5.6
#27    3       3 Compliance GW-2   11.0
#28    3       4 Compliance GW-2    9.8
#29    4       1 Compliance GW-2    6.3
#30    4       2 Compliance GW-2   10.4
#31    4       3 Compliance GW-2    7.5
#32    4       4 Compliance GW-2    9.7



# Step 1:  Summary Statistics and Time Series Plots
#--------------------------------------------------


# Summary Statistics by Well
#---------------------------

summaryStats(VC.ppb ~ Well, data = EPA.09.Ex.22.1.VC.df, 
	combine.groups = TRUE, digits = 1, stats.in.rows = TRUE)
#       GW-1 GW-2 Combined
#N      16   16   32      
#Mean    8.3  9    8.6    
#SD      1.9  2.9  2.4    
#Median  8.6  9.2  8.7    
#Min     4.9  3    3      
#Max    11.9 13.8 13.8 


# Summary Statistics by Period and Well
#--------------------------------------
VC     <- EPA.09.Ex.22.1.VC.df$VC.ppb
Period <- EPA.09.Ex.22.1.VC.df$Period
Well   <- EPA.09.Ex.22.1.VC.df$Well

summaryStats(VC, group = paste(Well, Period, sep = "."), 
	combine.groups = FALSE, stats.in.rows = TRUE, digits = 1)
#       GW-1.Background GW-1.Compliance GW-2.Background GW-2.Compliance
#N       4              12               4              12             
#Mean    8.9             8.1             7.4             9.5           
#SD      2.4             1.8             3.9             2.5           
#Median  8.8             8.6             7.4             9.8           
#Min     6.3             4.9             3               5.6           
#Max    11.9            11.2            12              13.8


# Group Shapiro-Wilk Test Shows p-values for Individual Wells
#------------------------------------------------------------

gofGroupTest(VC.ppb ~ Well, data = EPA.09.Ex.22.1.VC.df)
#
#Results of Group Goodness-of-Fit Test
#-------------------------------------
#
#Test Method:                     Wilk-Shapiro GOF (Normal Scores)
#
#Hypothesized Distribution:       Normal
#
#Data:                            VC.ppb
#
#Grouping Variable:               Well
#
#Data Source:                     EPA.09.Ex.22.1.VC.df
#
#Number of Groups:                2
#
#Sample Sizes:                    GW-1 = 16
#                                 GW-2 = 16
#
#Test Statistic:                  z (G) = 2.901581
#
#P-values for
#Individual Tests:                GW-1 = 0.9686679
#                                 GW-2 = 0.9875157
#
#P-value for
#Group Test:                      0.9981436
#
#Alternative Hypothesis:          At least one group
#                                 does not come from a
#                                 Normal Distribution.



# Time Series Plot
#-----------------
Year   <- EPA.09.Ex.22.1.VC.df$Year
Qtr    <- EPA.09.Ex.22.1.VC.df$Quarter

n1 <- sum(Well == "GW-1")
n2 <- sum(Well == "GW-2")
windows()
plot(1:n1, VC[Well == "GW-1"], xaxt = "n", 
	type = "b", pch = 15, lwd = 2, 
	ylim = c(0, 20),
	xlab = "Time (Years)",
	ylab = "Vinyl Chloride Concentration (ppb)",
	main = "Figure 22-1. Vinyl Chloride Time Series Plot")
	
axis(1, at = (1:n1)[Qtr == 1 & Well == "GW-1"], 
	labels = Year[Qtr == "1" & Well == "GW-1"])
points(1:n2, VC[Well == "GW-2"], type = "b", pch = 3, lty = 3, lwd = 2)
abline(h = 5, lty = 2, lwd = 2)
legend(1, 20, c("GW-1", "GW-2", "MCL = 5 ppb"), 
	pch = c(15, 3, 1), lty = c(1, 3, 2), 
	pt.cex = c(1, 1, 0), lwd = 2)



# Step 2:  Determine Type I Error level necessary in order to achieve
#          80% power to detect an increase in mean VC of 2 times
#          the 5 ppb GWPS based on n=4 observations.
#          Then compute one-sided lower confidence limits
#          for each well using the Year 2 data based on 
#          this Type I error level. 
#--------------------------------------------------------------------

R   <- 2
mu0 <- 5


# Approach 1:  
# Based on Guidance Document Method That *Does Not* Require
# an Estimate of SD or CV.  This method assumes:
# 1) sigma = mean under null hypothesis
# 2) delta/sigma = R - 1, where R is the ratio of the alternative mean
#    to the null mean
#------------------------------------------------------------------

required.alpha <- tTestAlpha(n.or.n1 = 4, 
	delta.over.sigma = R-1, power = 0.8, 
	sample.type = "one.sample", 
	alternative = "greater")
required.alpha
#[1] 0.1633382


# Lower one-sided confidence limits for Well GW-1:

enorm(VC[Well == "GW-1" & Year == "2"], ci = T, ci.type = "lower", 
	conf.level = 1 - required.alpha)$interval$limits["LCL"]
#     LCL 
#6.693358


# Lower one-sided confidence limits for Well GW-2:

enorm(VC[Well == "GW-2" & Year == "2"], ci = T, ci.type = "lower", 
	conf.level = 1 - required.alpha)$interval$limits["LCL"]
#     LCL 
#8.469282



# Approach 2:  
# Based on standard method that requires an estimate of SD or CV. 
#----------------------------------------------------------------

# Estimate standard deviation of VC based on Background data.
# Pool estimate of SDs based on both wells (GW-1 and GW-2)

VC.fit <- lm(VC.ppb ~ Well, data = EPA.09.Ex.22.1.VC.df, 
	subset = Period == "Background")
SD.Background <- summary(VC.fit)$sigma
SD.Background
#[1] 3.200976

#NOTE:  Note that the estimated SD of 3.2 ppb is much smaller than the 
#       assumed value of 5 ppb based on Approach 1 above.

tTestAlpha(n.or.n1 = 4, 
	delta.over.sigma = mu0 * (R - 1) / SD.Background, 
	power = 0.8, 
	sample.type = "one.sample", 
	alternative = "greater", approx = F)
#[1] 0.05767842

#NOTE:  Note the required Type I Error level based on Approach 2
#       is much smaller than that based on Approach 1.



# Step 3:  Compute one-sided lower confidence limits for each well
#          using the Year 2 data based on Type I Error Level = 0.01. 
#-------------------------------------------------------------------

# Use alpha = 0.01

# Lower one-sided confidence limits for Well GW-1:

enorm(VC[Well == "GW-1" & Year == "2"], ci = T, ci.type = "lower", 
	conf.level = 0.99)$interval$limits["LCL"]
#     LCL 
#2.926725


# Lower one-sided confidence limits for Well GW-2:

enorm(VC[Well == "GW-2" & Year == "2"], ci = T, ci.type = "lower", 
	conf.level = 0.99)$interval$limits["LCL"]
#    LCL 
#4.34498



# Step 4:  Determine Type I Error level necessary in order to achieve
#          80% power to detect an increase in mean VC of 2 times
#          the 5 ppb GWPS based on n=8 observations.
#          Then compute one-sided lower confidence limits
#          for each well using the Year 1 and Year 2 data based on 
#          this Type I error level. 
#--------------------------------------------------------------------

required.alpha <- tTestAlpha(n.or.n1 = 8, 
	delta.over.sigma = R-1, power = 0.8, 
	sample.type = "one.sample", 
	alternative = "greater")
required.alpha
#[1] 0.0458989

# Lower one-sided confidence limits for Well GW-1:

enorm(VC[Well == "GW-1" & Year %in% c("1", "2")], ci = T, ci.type = "lower", 
	conf.level = 1 - required.alpha)$interval$limits["LCL"]
#     LCL 
#6.963887


# Lower one-sided confidence limits for Well GW-2:

enorm(VC[Well == "GW-2" & Year %in% c("1", "2")], ci = T, ci.type = "lower", 
	conf.level = 1 - required.alpha)$interval$limits["LCL"]
#     LCL 
#6.403578 



# Step 5:  Compute one-sided lower confidence limits
#          for each well using the Year 2 and Year 3 data based on 
#          Type I error level from Step 4. 
#--------------------------------------------------------------------

# Lower one-sided confidence limits for Well GW-1:

enorm(VC[Well == "GW-1" & Year %in% c("2", "3")], ci = T, ci.type = "lower", 
	conf.level = 1 - required.alpha)$interval$limits["LCL"]
#     LCL 
#6.227279


# Lower one-sided confidence limits for Well GW-2:

enorm(VC[Well == "GW-2" & Year %in% c("2", "3")], ci = T, ci.type = "lower", 
	conf.level = 1 - required.alpha)$interval$limits["LCL"]
#     LCL 
#8.078256  



# Step 6:  Determine Type I Error level necessary in order to achieve
#          80% power to detect an increase in mean VC of 2 times
#          the 5 ppb GWPS based on n=12 observations.
#          Then compute one-sided lower confidence limits
#          for each well using the Years 2-4 data based on 
#          this Type I error level. 
#--------------------------------------------------------------------

required.alpha <- tTestAlpha(n.or.n1 = 12, 
	delta.over.sigma = R-1, power = 0.8, 
	sample.type = "one.sample", 
	alternative = "greater")
required.alpha
#[1] 0.01311214


# Lower one-sided confidence limits for Well GW-1:

enorm(VC[Well == "GW-1" & Year != "1"], ci = T, ci.type = "lower", 
	conf.level = 1 - required.alpha)$interval$limits["LCL"]
#     LCL 
#6.798091


# Lower one-sided confidence limits for Well GW-2:

enorm(VC[Well == "GW-2" & Year != "1"], ci = T, ci.type = "lower", 
	conf.level = 1 - required.alpha)$interval$limits["LCL"]
#     LCL 
#7.609904 



rm(Year, Qtr, Period, Well, VC, n1, n2, 
	R, mu0, required.alpha, VC.fit, SD.Background)



#########################################################################################


# Example 22-2, pp. 22-11 to 22-17
#---------------------------------

EPA.09.Ex.22.2.Specific.Conductance.df
#    Well       Date Specific.Conductance.umho
#1  GW-12 1987-10-16                      2100
#2  GW-12 1988-01-28                      2550
#   ...
#42 GW-13 1992-12-01                       530
#43 GW-13 1993-03-24                       625

Well <- EPA.09.Ex.22.2.Specific.Conductance.df$Well
Date <- EPA.09.Ex.22.2.Specific.Conductance.df$Date
SC   <- EPA.09.Ex.22.2.Specific.Conductance.df$Specific.Conductance.umho


#Step 1:  Time series plot for Well GW-12
#---------------------------------------

index <- Well == "GW-12"

windows()
plot(Date[index], SC[index], ylim = c(0, 3000), 
	xlab = "Sampling Date", 
	ylab = "Specific Conductance (umho)",
	main = "Figure 22-3. Time Series Plot of \nSpecific Conductance Measurements at GW-12")


#Step 2:  Fit regression line
#----------------------------

SC.fit <- lm(Specific.Conductance.umho ~ Date, 
	data = EPA.09.Ex.22.2.Specific.Conductance.df, 
	subset = Well == "GW-12")

summary(SC.fit)$coef
#               Estimate   Std. Error  t value     Pr(>|t|)
#(Intercept) 9526.116147 1042.1592274  9.14075 3.493624e-08
#Date          -1.087567    0.1384346 -7.85618 3.170120e-07

# NOTE:  These estimates of slope and intercept differ from those in the 
#        Guidance Document because our time variable is in terms of 
#        days since Jan 1, 1970, where as the units of time in the 
#        Guidance Document is Years.


windows()
plot(Date[index], SC[index], ylim = c(0, 3000), 
	xlab = "Sampling Date", 
	ylab = "Specific Conductance (umho)",
	main = "Figure 22-4. Regression of Specific Conductance vs.\nSampling Date at GW-12")
abline(SC.fit)


# Step 3:  Examine residuals
#---------------------------

windows()
qqPlot(SC.fit$residuals, xlab = "Z-score", add.line = T, 
	ylab = "Specific Conductance Residuals",
	main = "Figure 22-5. Probability Plot of \nSpecific Conductance Residuals at GW-12")

gofTest(SC.fit$residuals)
#
#Results of Goodness-of-Fit Test
#-------------------------------
#
#Test Method:                     Shapiro-Wilk GOF
#
#Hypothesized Distribution:       Normal
#
#Estimated Parameter(s):          mean = 3.538836e-16
#                                 sd   = 3.790014e+02
#
#Estimation Method:               mvue
#
#Data:                            SC.fit$residuals
#
#Sample Size:                     20
#
#Test Statistic:                  W = 0.960157
#
#Test Statistic Parameter:        n = 20
#
#P-value:                         0.5469918
#
#Alternative Hypothesis:          True cdf does not equal the
#                                 Normal Distribution.

windows()
plot(Date[index], SC.fit$residuals, xlab = "Sampling Date", 
	ylab = "Specific Conductance Residuals",
	main = "Figure 22-6. Plot of Specific Conductance Residuals vs. \nSampling Date at GW-12")
abline(h = 0, lty = 2)


# Alternatively, you can use used the built-in plotting function for linear models:
#windows()
#plot(SC.fit, ask = FALSE)



# Step 4:  Plot Specific Conductance at GW-12 vs. Date, along with
#          regression line and two-sided lower and upper 90% confidence limits
#          Note that the upper confidence limit considered by itself is 
#          equivalent to a one-sided upper 95% confidence limit.
#-----------------------------------------------------------------------------

Dates.for.prediction <- seq(min(Date[index]), max(Date[index]), length = 100)
CLs <- predict(SC.fit, newdata = data.frame(Date = Dates.for.prediction),
	interval = "confidence", level = 0.9)

windows()
plot(Date[index], SC[index], ylim = c(0, 3000), 
	xlab = "Sampling Date", 
	ylab = "Specific Conductance (umho)",
	main = "Figure 22-7. Two-sided 90% Confidence Bounds \nAround Trend Line at GW-12")
abline(SC.fit, lwd = 2)
lines(Dates.for.prediction, CLs[, "lwr"], lty = 3, lwd = 2, col = 4)
lines(Dates.for.prediction, CLs[, "upr"], lty = 3, lwd = 2, col = 4)
abline(h = 1000, lty = 2, lwd = 2, col = 2)

legend(as.numeric(as.Date("1988-01-01")), 500, 
	c("Regression Line", "90% Conf Limits", "Clean-Up Std"),
	lty = c(1, 3, 2), lwd = 2, col = c(1, 4, 2))


# Step 5:  Determine when the One-sided upper 95% confidence limit
#          drops below the clean-up standard of 1000 umho.

min(Dates.for.prediction[CLs[, "upr"] < 1000])
#[1] "1991-12-18".




# Step 6: Time series plot for Well GW-13
#         Also, remove first 4 observations
#------------------------------------------

index <- Well == "GW-13"

windows()
plot(Date[index], SC[index], type = "b", 
	ylim = c(0, 2500), 
	xlab = "Sampling Date", 
	ylab = "Specific Conductance (umho)",
	main = "Figure 22-8. Time Series Plot of \nSpecific Conductance Measurements at GW-13")
 
new.index <- index & (Date > Date[index][4])



# Step 7: Test normality of GW-13 data
#-------------------------------------

# Using all the data

gof.obj <- gofTest(SC[index])
gof.obj[c("statistic", "p.value")]
#$statistic
#        W 
#0.5809712 
#
#$p.value
#[1] 5.617319e-07

windows()
plot(gof.obj, plot.type = "Q-Q Plot",
	ylab = "Specific Conductance (umho)", 
	main = "Figure 22-9. Probability Plot at GW-13\nUsing Entire Sampling Record")


# Omitting the first 4 observations
#----------------------------------
gof.obj <- gofTest(SC[new.index])
gof.obj[c("statistic", "p.value")]
#$statistic
#        W 
#0.9542246 
#
#$p.value
#[1] 0.4647555

windows()
plot(gof.obj, plot.type = "Q-Q Plot",
	ylab = "Specific Conductance (umho)", 
	main = "Figure 22-10. Probability Plot at GW-13\nExcluding First Four Measures")


# Step 8:  Compute One-sided upper 95% confidence bounds for 
#          mean Specific Conductance at Well GW-13 using the first 8 
#          observations, and then omitting the first 4 observations
#--------------------------------------------------------------------

# Using the first 8 observations
enorm(SC[index][1:8], ci = T, ci.type = "upper")$interval$limits["UCL"]
#     UCL 
#1327.572 


# Using observations 5 to 8
enorm(SC[index][5:8], ci = T, ci.type = "upper")$interval$limits["UCL"]
#     UCL 
#515.3008



# Step 9:  Compute One-sided upper 95% confidence bounds for 
#          mean Specific Conductance at Well GW-13 using 
#          observations 5 to 12, and then using all observations 
#          except the first 4.
#--------------------------------------------------------------------


# Using observations 5 to 12
enorm(SC[index][5:12], ci = T, ci.type = "upper")$interval$limits["UCL"]
#   UCL 
#524.53


# Using all observations except the first 4
enorm(SC[new.index], ci = T, ci.type = "upper")$interval$limits["UCL"]
#     UCL 
#532.7065 



# Step 10:  Compute sample size necessary to detect a decrease in 
#           mean conductance from 1000 umho to 500 umho
#           based on estimated SD using all the data at Well GW-13
#           except the first 4 observations.
#-----------------------------------------------------------------

SD <- sd(SC[new.index])
SD
#[1] 74.27596

tTestN(delta.over.sigma = (500 - 1000) / SD, 
	alpha = 0.05, power = 0.95, 
	sample.type = "one.sample", alternative = "less", round.up = T)
#[1] 3

# NOTE:  The Guidance Document states n=6 is necessary using 
#        it's second approach.


rm(Well, Date, SC, index, SC.fit, Dates.for.prediction, 
	CLs, new.index, gof.obj, SD)



###############################################################################################


# Example 22-3, p. 22-20
#-----------------------


# Based on normal approximation to the binomial
#----------------------------------------------
propTestN(p.or.p1 = 0.3, p0.or.p2 = 0.1, alpha = 0.05, power = 0.9, 
	sample.type = "one.sample", alternative = "greater",   
	approx = T, round.up = T)
#[1] 30


# Based on exact binomial distribution
#-------------------------------------
propTestN(p.or.p1 = 0.3, p0.or.p2 = 0.1, alpha = 0.05, power = 0.9, 
	sample.type = "one.sample", alternative = "greater",   
	approx = F, round.up = T)
#$n
#[1] 34
#
#$power
#[1] 0.9214717
#
#$alpha
#[1] 0.04814433
#
#$q.critical.upper
#[1] 6


propTestN(p.or.p1 = 0.3, p0.or.p2 = 0.1, alpha = 0.07, power = 0.9,
	sample.type = "one.sample", alternative = "greater",   
	approx = F, round.up = T)
#$n
#[1] 29
#
#$power
#[1] 0.9068196
#
#$alpha
#[1] 0.06371744
#
#$q.critical.upper
#[1] 5


propTestN(p.or.p1 = 0.3, p0.or.p2 = 0.1, alpha = 0.1, power = 0.9, 
	sample.type = "one.sample", alternative = "greater",   
	approx = F, round.up = T)
#$n
#[1] 25
#
#$power
#[1] 0.909528
#
#$alpha
#[1] 0.09799362
#
#$q.critical.upper
#[1] 4


# So using the exact Binomial test, you would need a minimum sample size of 34, NOT 30,
# to detect an increase from 10% to 30% with 90% power and a Type I error of
# no greater than 5% (4.8%) in this case.
# Allowing a Type I Error of 6.3% would let us use a sample size of n=29. 


############################################################################################


# Example 22-4, pp. 22-21 to 22-22
#---------------------------------

# Step 1: Compute sample size necessary to detect a change in proprtion of
#         exceedances from 0.05 to (0.25 * 0.05) = 0.0125 using
#         a Type I error of 20% or 10% and a power of 80%
#--------------------------------------------------------------------------


k <- 0.25
p0 <- 0.05

# Based on normal approximation to the binomial
#----------------------------------------------

propTestN(p.or.p1 = k * p0, p0.or.p2 = p0, 
	alpha = 0.2, power = 0.8, 
	sample.type = "one.sample", alternative = "less",   
	approx = T, round.up = T)
#[1] 55
#Warning message:
#In propTestN(p.or.p1 = k * p0, p0.or.p2 = p0, alpha = 0.2,  :
#  The computed sample size 'n' is too small, relative to the given value of 'p', for the normal approximation to work well for the following element indices:
#         1 
    

propTestN(p.or.p1 = k * p0, p0.or.p2 = p0, 
	alpha = 0.1, power = 0.8, 
	sample.type = "one.sample", alternative = "less",   
	approx = T, round.up = T)
#[1] 99
#Warning message:
#In propTestN(p.or.p1 = k * p0, p0.or.p2 = p0, alpha = 0.1,  :
#  The computed sample size 'n' is too small, relative to the given value of 'p', for the normal approximation to work well for the following element indices:
#         1 




# Step 2a:  Compute sample size necessary to detect a change in proprtion of
#           exceedances from 0.1 to (0.25 * 0.1) = 0.025 using 
#           a Type I error of 10% and a power of 80%
#----------------------------------------------------------------------------

p0 <- 0.1

propTestN(p.or.p1 = k * p0, p0.or.p2 = p0, 
	alpha = 0.1, power = 0.8, 
	sample.type = "one.sample", alternative = "less",   
	approx = T, round.up = T)
#[1] 48
#Warning message:
#In propTestN(p.or.p1 = k * p0, p0.or.p2 = p0, alpha = 0.1, power = 0.8,  :
#  The computed sample size 'n' is too small, relative to the given value of 'p', for the normal approximation to work well for the following element indices:
#         1 


# Step 2b:  Compute sample size necessary to detect a change in proprtion of
#           exceedances from 0.05 to (0.1 * 0.05) = 0.005 using 
#           a Type I error of 10% and a power of 80%
#----------------------------------------------------------------------------

k <- 0.1
p0 <- 0.05

propTestN(p.or.p1 = k * p0, p0.or.p2 = p0, 
	alpha = 0.1, power = 0.8, 
	sample.type = "one.sample", alternative = "less",   
	approx = T, round.up = T)
#[1] 57
#Warning message:
#In propTestN(p.or.p1 = k * p0, p0.or.p2 = p0, alpha = 0.1, power = 0.8,  :
#  The computed sample size 'n' is too small, relative to the given value of 'p', for the normal approximation to work well for the following element indices:
#         1 


# Step 2c:  Compute sample size necessary to detect a change in proprtion of
#           exceedances from 0.05 to (0.1 * 0.05) = 0.005 using 
#           a Type I error of 10% and a power of 60%
#----------------------------------------------------------------------------

k  <- 0.25
p0 <- 0.05

propTestN(p.or.p1 = k * p0, p0.or.p2 = p0, 
	alpha = 0.1, power = 0.6, 
	sample.type = "one.sample", alternative = "less",   
	approx = T, round.up = T)
#[1] 68
#Warning message:
#In propTestN(p.or.p1 = k * p0, p0.or.p2 = p0, alpha = 0.1, power = 0.6,  :
#  The computed sample size 'n' is too small, relative to the given value of 'p', for the normal approximation to work well for the following element indices:
#         1 


rm(k, p0)


