#####
# File:     Script - EPA09, Chapter 21 Examples.R
#
# Purpose:  Reproduce Examples in Chapter 21 of the EPA Guidance Document
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

# Example 21-1, pp. 21-4 to 21-5
#-------------------------------

# Raw data
#---------

EPA.09.Ex.21.1.aldicarb.df
#   Month   Well Aldicarb.ppb
#1      1 Well.1         19.9
#2      2 Well.1         29.6
#3      3 Well.1         18.7
#4      4 Well.1         24.2
#5      1 Well.2         23.7
#6      2 Well.2         21.9
#7      3 Well.2         26.9
#8      4 Well.2         26.1
#9      1 Well.3          5.6
#10     2 Well.3          3.3
#11     3 Well.3          2.3
#12     4 Well.3          6.9

longToWide(EPA.09.Ex.21.1.aldicarb.df, "Aldicarb.ppb", 
	"Month", "Well", paste.row.name = TRUE)
#        Well.1 Well.2 Well.3
#Month.1   19.9   23.7    5.6
#Month.2   29.6   21.9    3.3
#Month.3   18.7   26.9    2.3
#Month.4   24.2   26.1    6.9



# Step 1:  Summary Statistics and Shapiro-Wilk Tests 
#---------------------------------------------------


# Summary Statistics by Well
#---------------------------

summaryStats(Aldicarb.ppb ~ Well, 
	data = EPA.09.Ex.21.1.aldicarb.df, 
	combine.groups = TRUE, digits = 1, stats.in.rows = TRUE)
#       Well.1 Well.2 Well.3 Combined
#N       4      4      4     12      
#Mean   23.1   24.6    4.5   17.4    
#SD      4.9    2.3    2.1   10      
#Median 22     24.9    4.4   20.9    
#Min    18.7   21.9    2.3    2.3    
#Max    29.6   26.9    6.9   29.6 


summaryFull(Aldicarb.ppb ~ Well, data = EPA.09.Ex.21.1.aldicarb.df,
	stats = c("n", "mean", "sd", "skew"))
#                   Well.1  Well.2  Well.3 
#N                   4       4       4     
#Mean               23.1    24.65    4.525 
#Skew                0.8765 -0.4045  0.1283
#Standard Deviation  4.935   2.283   2.101



# Group Shapiro-Wilk Test Shows p-values for Individual Wells
#------------------------------------------------------------

gofGroupTest(Aldicarb.ppb ~ Well, data = EPA.09.Ex.21.1.aldicarb.df)
#
#Results of Group Goodness-of-Fit Test
#-------------------------------------
#
#Test Method:                     Wilk-Shapiro GOF (Normal Scores)
#
#Hypothesized Distribution:       Normal
#
#Data:                            Aldicarb.ppb
#
#Grouping Variable:               Well
#
#Data Source:                     EPA.09.Ex.21.1.aldicarb.df
#
#Number of Groups:                3
#
#Sample Sizes:                    Well.1 = 4
#                                 Well.2 = 4
#                                 Well.3 = 4
#
#Test Statistic:                  z (G) = 0.618909
#
#P-values for
#Individual Tests:                Well.1 = 0.5469869
#                                 Well.2 = 0.6617847
#                                 Well.3 = 0.7042248
#
#P-value for
#Group Test:                      0.7320118
#
#Alternative Hypothesis:          At least one group
#                                 does not come from a
#                                 Normal Distribution.



# Steps 2-4 Compute 95% Lower Confidence Bounds for the mean 
# and Compare to the standard of 7 ppb 
#-----------------------------------------------------------

# Example for Well 1
#-------------------

Month    <- EPA.09.Ex.21.1.aldicarb.df$Month
Well     <- EPA.09.Ex.21.1.aldicarb.df$Well
Aldicarb <- EPA.09.Ex.21.1.aldicarb.df$Aldicarb.ppb

enorm(Aldicarb[Well == "Well.1"], ci = T, ci.type = "lower")
#
#Results of Distribution Parameter Estimation
#--------------------------------------------
#
#Assumed Distribution:            Normal
#
#Estimated Parameter(s):          mean = 23.10000
#                                 sd   =  4.93491
#
#Estimation Method:               mvue
#
#Data:                            Aldicarb[Well == "Well.1"]
#
#Sample Size:                     4
#
#Confidence Interval for:         mean
#
#Confidence Interval Method:      Exact
#
#Confidence Interval Type:        lower
#
#Confidence Level:                95%
#
#Confidence Interval:             LCL = 17.29318
#                                 UCL =      Inf


# Example for all 3 Wells at Once
#--------------------------------

LCLs <- sapply(split(Aldicarb, Well), 
	function(x) enorm(x, ci = T, ci.type = "lower")$interval$limits["LCL"])
round(LCLs, 2)
#Well.1.LCL Well.2.LCL Well.3.LCL 
#     17.29      21.96       2.05 

LCLs > 7
#Well.1.LCL Well.2.LCL Well.3.LCL 
#      TRUE       TRUE      FALSE


rm(Month, Well, Aldicarb, LCLs)


##########################################################################################


# Example 21-2, pp. 21-7 to 21-8
#-------------------------------

# Raw Data
#---------

EPA.09.Ex.21.2.benzene.df
#  Month Benzene.ppb.orig Benzene.ppb Censored
#1     1              0.5         0.5    FALSE
#2     2              0.5         0.5    FALSE
#3     3              1.6         1.6    FALSE
#4     4              1.8         1.8    FALSE
#5     5              1.1         1.1    FALSE
#6     6             16.1        16.1    FALSE
#7     7              1.6         1.6    FALSE
#8     8             <0.5         0.5     TRUE


# Step 1:  Summary Statistics and Shapiro-Wilk Tests 
#---------------------------------------------------

# Impute non-detect of 0.5 ppb as 0.25 ppb
Benzene      <- EPA.09.Ex.21.2.benzene.df$Benzene.ppb
Censored     <- EPA.09.Ex.21.2.benzene.df$Censored
Benzene[Censored] <- 0.25


# Summary Statistics for raw data
#--------------------------------
summaryStats(Benzene)
#        N Mean  SD Median Min  Max
#Benzene 8  2.9 5.4    1.4 0.2 16.1
 
summaryFull(Benzene, stats = c("n", "mean", "sd", "skew"))
#                    Benzene
#Sample Size:        8      
#Mean:               2.931  
#Skew:               2.76   
#Standard Deviation: 5.353


# Shapiro-Wilk Test on raw data
#------------------------------
gofTest(Benzene)
#
#Results of Goodness-of-Fit Test
#-------------------------------
#
#Test Method:                     Shapiro-Wilk GOF
#
#Hypothesized Distribution:       Normal
#
#Estimated Parameter(s):          mean = 2.931250
#                                 sd   = 5.353299
#
#Estimation Method:               mvue
#
#Data:                            Benzene
#
#Sample Size:                     8
#
#Test Statistic:                  W = 0.5208992
#
#Test Statistic Parameter:        n = 8
#
#P-value:                         1.841521e-05
#
#Alternative Hypothesis:          True cdf does not equal the
#                                 Normal Distribution.



# Summary Statistics for log-transformed data
#--------------------------------------------
log.Benzene <- log(Benzene)

summaryStats(log.Benzene)
#            N Mean  SD Median  Min Max
#log.Benzene 8  0.2 1.3    0.3 -1.4 2.8
 
summaryFull(log.Benzene, stats = c("n", "mean", "sd", "skew"))
#                    log.Benzene
#Sample Size:        8          
#Mean:               0.2037     
#Skew:               1.1215     
#Standard Deviation: 1.2575


# Shapiro-Wilk Test for lognormal distribution
#---------------------------------------------
gofTest(Benzene, dist = "lnorm")
#
#Results of Goodness-of-Fit Test
#-------------------------------
#
#Test Method:                     Shapiro-Wilk GOF
#
#Hypothesized Distribution:       Lognormal
#
#Estimated Parameter(s):          meanlog = 0.2036668
#                                 sdlog   = 1.2574974
#
#Estimation Method:               mvue
#
#Data:                            Benzene
#
#Sample Size:                     8
#
#Test Statistic:                  W = 0.896502
#
#Test Statistic Parameter:        n = 8
#
#P-value:                         0.2686297
#
#Alternative Hypothesis:          True cdf does not equal the
#                                 Lognormal Distribution.


# Steps 2-4 Compute 95% Upper Confidence Bound for 
# the Geometric Mean and compare to MCL of 5 ppb 
#-------------------------------------------------

exp(elnorm(Benzene, ci = T, ci.type = "upper")$interval$limits["UCL"])
#     UCL 
#2.846193 

########################################################################################

# Example 21-3, pp. 21-10 to 21-11
#---------------------------------

# Using the same data as in Example 21-2, 
# compute 95% Upper Confidence Bound for Mean

elnormAlt(Benzene, ci = T, ci.type = "upper")
#
#Results of Distribution Parameter Estimation
#--------------------------------------------
#
#Assumed Distribution:            Lognormal
#
#Estimated Parameter(s):          mean = 2.338968
#                                 cv   = 1.331584
#
#Estimation Method:               mvue
#
#Data:                            Benzene
#
#Sample Size:                     8
#
#Confidence Interval for:         mean
#
#Confidence Interval Method:      Land
#
#Confidence Interval Type:        upper
#
#Confidence Level:                95%
#
#Confidence Interval:             LCL =  0.00000
#                                 UCL = 18.86139


rm(Benzene, Censored, log.Benzene)


##########################################################################################


# Example 21-4, pp. 21-13 to 21-14
#---------------------------------

# Raw Data - same as in Example 21-1
#-----------------------------------

EPA.09.Ex.21.1.aldicarb.df
#   Month   Well Aldicarb.ppb
#1      1 Well.1         19.9
#2      2 Well.1         29.6
#3      3 Well.1         18.7
#4      4 Well.1         24.2
#5      1 Well.2         23.7
#6      2 Well.2         21.9
#7      3 Well.2         26.9
#8      4 Well.2         26.1
#9      1 Well.3          5.6
#10     2 Well.3          3.3
#11     3 Well.3          2.3
#12     4 Well.3          6.9

longToWide(EPA.09.Ex.21.1.aldicarb.df, "Aldicarb.ppb", 
	"Month", "Well", paste.row.name = TRUE)
#        Well.1 Well.2 Well.3
#Month.1   19.9   23.7    5.6
#Month.2   29.6   21.9    3.3
#Month.3   18.7   26.9    2.3
#Month.4   24.2   26.1    6.9



# Compute 99% Lower Confidence Limits for the 95th Percentile
# of Aldicarb for Compliance/Assessment Monitoring
#------------------------------------------------------------

# Example for Well 1
#-------------------

Month    <- EPA.09.Ex.21.1.aldicarb.df$Month
Well     <- EPA.09.Ex.21.1.aldicarb.df$Well
Aldicarb <- EPA.09.Ex.21.1.aldicarb.df$Aldicarb.ppb

eqnorm(Aldicarb[Well == "Well.1"], p = 0.95, method = "qmle", ci = T, 
	ci.type = "lower", conf.level = 0.99)
#
#Results of Distribution Parameter Estimation
#--------------------------------------------
#
#Assumed Distribution:            Normal
#
#Estimated Parameter(s):          mean = 23.10000
#                                 sd   =  4.93491
#
#Estimation Method:               mvue
#
#Estimated Quantile(s):           95'th %ile = 31.21720
#
#Quantile Estimation Method:      qmle
#
#Data:                            Aldicarb[Well == "Well.1"]
#
#Sample Size:                     4
#
#Confidence Interval for:         95'th %ile
#
#Confidence Interval Method:      Exact
#
#Confidence Interval Type:        lower
#
#Confidence Level:                99%
#
#Confidence Interval:             LCL = 25.28550
#                                 UCL =      Inf
#


# Example for all 3 Wells at Once
#--------------------------------

LCLs <- sapply(split(Aldicarb, Well), 
	function(x) eqnorm(x, p = 0.95, method = "qmle", ci = T, 
	ci.type = "lower", conf.level = 0.99)$interval$limits["LCL"])
round(LCLs, 2)
#Well.1.LCL Well.2.LCL Well.3.LCL 
#     25.29      25.66       5.46 

LCLs > 30
#Well.1.LCL Well.2.LCL Well.3.LCL 
#     FALSE      FALSE      FALSE



# Compute 99% Upper Confidence Limits for the 95th Percentile
# of Aldicarb for Corrective Action Monitoring
#------------------------------------------------------------

# Example for Well 1
#-------------------

eqnorm(Aldicarb[Well == "Well.1"], p = 0.95, method = "qmle", ci = T, 
	ci.type = "upper", conf.level = 0.99)
#
#Results of Distribution Parameter Estimation
#--------------------------------------------
#
#Assumed Distribution:            Normal
#
#Estimated Parameter(s):          mean = 23.10000
#                                 sd   =  4.93491
#
#Estimation Method:               mvue
#
#Estimated Quantile(s):           95'th %ile = 31.21720
#
#Quantile Estimation Method:      qmle
#
#Data:                            Aldicarb[Well == "Well.1"]
#
#Sample Size:                     4
#
#Confidence Interval for:         95'th %ile
#
#Confidence Interval Method:      Exact
#
#Confidence Interval Type:        upper
#
#Confidence Level:                99%
#
#Confidence Interval:             LCL =     -Inf
#                                 UCL = 67.92601


# Example for all 3 Wells at Once
#--------------------------------

UCLs <- sapply(split(Aldicarb, Well), 
	function(x) eqnorm(x, p = 0.95, method = "qmle", ci = T, 
	ci.type = "upper", conf.level = 0.99)$interval$limits["UCL"])
round(UCLs, 2)
#Well.1.UCL Well.2.UCL Well.3.UCL 
#     67.93      45.38      23.61 
 
UCLs > 30
#Well.1.UCL Well.2.UCL Well.3.UCL 
#      TRUE       TRUE      FALSE 


rm(Month, Well, Aldicarb, LCLs, UCLs)


#############################################################################################


# Example 21-5, pp. 21-18 to 21-20
#---------------------------------

EPA.09.Ex.21.5.beryllium.df
#   Year Quarter Beryllium.ppb
#1  2002       1          3.17
#2  2002       2          2.32
#3  2002       3          7.37
#4  2002       4          4.44
#5  2003       1          9.50
#6  2003       2         21.36
#7  2003       3          5.15
#8  2003       4         15.70
#9  2004       1          5.58
#10 2004       2          3.39
#11 2004       3          8.44
#12 2004       4         10.25
#13 2005       1          3.65
#14 2005       2          6.15
#15 2005       3          6.94
#16 2005       4          3.74

Be <- EPA.09.Ex.21.5.beryllium.df$Beryllium.ppb


# Steps 1-3:  Nonparametric lower 99% confidence limit for the median
#--------------------------------------------------------------------

eqnpar(Be, p = 0.5, ci = T, ci.type = "lower", approx.conf.level = 0.99)
#
#Results of Distribution Parameter Estimation
#--------------------------------------------
#
#Assumed Distribution:            None
#
#Estimated Quantile(s):           Median = 5.865
#
#Quantile Estimation Method:      Nonparametric
#
#Data:                            Be
#
#Sample Size:                     16
#
#Confidence Interval for:         50'th %ile
#
#Confidence Interval Method:      exact
#
#Confidence Interval Type:        lower
#
#Confidence Level:                99%
#
#Confidence Limit Rank(s):        4 
#
#Confidence Interval:             LCL = 3.65
#                                 UCL =  Inf

# Note, if you want to see the confidence level to 2 digits, do this:

eqnpar.obj <- eqnpar(Be, p = 0.5, ci = T, ci.type = "lower", approx.conf.level = 0.99)
print(eqnpar.obj, conf.cov.sig.digits=4)
#
#Results of Distribution Parameter Estimation
#--------------------------------------------
#
#Assumed Distribution:            None
#
#Estimated Quantile(s):           Median = 5.865
#
#Quantile Estimation Method:      Nonparametric
#
#Data:                            Be
#
#Sample Size:                     16
#
#Confidence Interval for:         50'th %ile
#
#Confidence Interval Method:      exact
#
#Confidence Interval Type:        lower
#
#Confidence Level:                98.94%
#
#Confidence Limit Rank(s):        4 
#
#Confidence Interval:             LCL = 3.65
#                                 UCL =  Inf


# Step 4:  Compute 95% Lower Confidence Bound for 
# the Geometric Mean and compare to Nonparametric LCL 
#----------------------------------------------------

exp(elnorm(Be, ci = T, ci.type = "lower", conf.level = 0.99)$interval$limits["LCL"])
#     LCL 
#4.129290

rm(Be, eqnpar.obj)

#######################################################################################


# Example 21-6, pp. 21-21 to 21-22
#---------------------------------

EPA.09.Ex.21.6.nitrate.df
#   Sampling.Date       Date Nitrate.mg.per.l.orig Nitrate.mg.per.l Censored
#1      7/28/1999 1999-07-28                  <5.0              5.0     TRUE
#2       9/3/1999 1999-09-03                  12.3             12.3    FALSE
#3     11/24/1999 1999-11-24                  <5.0              5.0     TRUE
#4       5/3/2000 2000-05-03                  <5.0              5.0     TRUE
#5      7/14/2000 2000-07-14                   8.1              8.1    FALSE
#6     10/31/2000 2000-10-31                  <5.0              5.0     TRUE
#7     12/14/2000 2000-12-14                    11             11.0    FALSE
#8      3/27/2001 2001-03-27                  35.1             35.1    FALSE
#9      6/13/2001 2001-06-13                  <5.0              5.0     TRUE
#10     9/16/2001 2001-09-16                  <5.0              5.0     TRUE
#11    11/26/2001 2001-11-26                   9.3              9.3    FALSE
#12      3/2/2002 2002-03-02                  10.3             10.3    FALSE



Nitrate <- EPA.09.Ex.21.6.nitrate.df$Nitrate.mg.per.l


# Steps 1-4:  Compute nonparametric 95% LCL for 95th percentile
#--------------------------------------------------------------


eqnpar(Nitrate, p = 0.95, ci = T, ci.type = "lower", approx.conf.level = 0.95)
#
#Results of Distribution Parameter Estimation
#--------------------------------------------
#
#Assumed Distribution:            None
#
#Estimated Quantile(s):           95'th %ile = 22.56
#
#Quantile Estimation Method:      Nonparametric
#
#Data:                            Nitrate
#
#Sample Size:                     12
#
#Confidence Interval for:         95'th %ile
#
#Confidence Interval Method:      exact
#
#Confidence Interval Type:        lower
#
#Confidence Level:                88%
#
#Confidence Limit Rank(s):        11 
#
#Confidence Interval:             LCL = 12.3
#                                 UCL =  Inf


#NOTE:  achieved confidence level is only 88%, so need to either increase
#       the value of approx.conf.level or else set lcl.rank to 10.

eqnpar(Nitrate, p = 0.95, ci = T, ci.type = "lower", lcl.rank = 10)
#
#Results of Distribution Parameter Estimation
#--------------------------------------------
#
#Assumed Distribution:            None
#
#Estimated Quantile(s):           95'th %ile = 22.56
#
#Quantile Estimation Method:      Nonparametric
#
#Data:                            Nitrate
#
#Sample Size:                     12
#
#Confidence Interval for:         95'th %ile
#
#Confidence Interval Method:      exact
#
#Confidence Interval Type:        lower
#
#Confidence Level:                98%
#
#Confidence Limit Rank(s):        10 
#
#Confidence Interval:             LCL =  11
#                                 UCL = Inf


# Compare LCL to acute risk standard of 10 mg/L
eqnpar.obj <- eqnpar(Nitrate, p = 0.95, ci = T, ci.type = "lower", lcl.rank = 10)
eqnpar.obj$interval$limits['LCL']
#LCL 
# 11 

 
eqnpar.obj$interval$limits['LCL'] < 10
#  LCL 
#FALSE 


# Step 5:  Compare nonparametric Upper 95% confidence limit for 95th percentile
#          to a corrective action limit of 10 ppm
#------------------------------------------------------------------------------

eqnpar(Nitrate, p = 0.95, ci = T, ci.type = "upper", approx.conf.level = 0.95)
#
#Results of Distribution Parameter Estimation
#--------------------------------------------
#
#Assumed Distribution:            None
#
#Estimated Quantile(s):           95'th %ile = 22.56
#
#Quantile Estimation Method:      Nonparametric
#
#Data:                            Nitrate
#
#Sample Size:                     12
#
#Confidence Interval for:         95'th %ile
#
#Confidence Interval Method:      exact
#
#Confidence Interval Type:        upper
#
#Confidence Level:                46%
#
#Confidence Limit Rank(s):        12 
#
#Confidence Interval:             LCL = -Inf
#                                 UCL = 35.1

#NOTE: Even using the maximum value as the UCL for the 95'th percentile, 
#      only a confidence level of 46% is achieved.


rm(Nitrate, eqnpar.obj)

#############################################################################################


# Example 21-7, pp. 21-26 to 21-29
#---------------------------------

EPA.09.Ex.21.7.TCE.df
#   Month TCE.ppb
#1      2    54.2
#2      4    44.3
#3      8    45.4
#4     11    38.3
#5     13    27.1
#6     16    30.2
#7     20    28.3
#8     23    17.6
#9     26    14.7
#10    30     4.1

Month <- EPA.09.Ex.21.7.TCE.df$Month
TCE   <- EPA.09.Ex.21.7.TCE.df$TCE.ppb

# Step 1 - Time series plot and regression line
#----------------------------------------------

windows()
plot(Month, TCE, xlim = c(0, 35), ylim = c(0, 60), 
	xlab = "Month Sampled", ylab = "TCE Conc (ppb)",
	main = "Figure 21-3. Time Series Plot and Regression Line of TCE Measurements")
TCE.fit <- lm(TCE ~ Month, data = EPA.09.Ex.21.7.TCE.df)
abline(TCE.fit)

summary(TCE.fit)
#
#Call:
#lm(formula = TCE ~ Month, data = EPA.09.Ex.21.7.TCE.df)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-7.0061 -2.1908  0.9452  2.2057  5.4124 
#
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  54.9405     2.4816   22.14 1.83e-08 ***
#Month        -1.6026     0.1402  -11.44 3.09e-06 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
#
#Residual standard error: 3.95 on 8 degrees of freedom
#Multiple R-squared: 0.9423,     Adjusted R-squared: 0.9351 
#F-statistic: 130.8 on 1 and 8 DF,  p-value: 3.094e-06 

# NOTE:  Estimated slope is -1.6026.  
#        Estimated standard error based on residuals is 3.95, so
#        Estimated variance based on residuals is 3.95^2 = 15.60



# Step 2 - Look at residuals from fit
#------------------------------------
windows()
plot(TCE.fit, which = 2, xlab = "Residual TCE (ppb)", 
	caption = "", 
	main = "Figure 21-4. Probability Plot of TCE Residuals")

gofTest(TCE.fit$residuals)
#
#Results of Goodness-of-Fit Test
#-------------------------------
#
#Test Method:                     Shapiro-Wilk GOF
#
#Hypothesized Distribution:       Normal
#
#Estimated Parameter(s):          mean = -1.332918e-16
#                                 sd   =  3.723714e+00
#
#Estimation Method:               mvue
#
#Data:                            TCE.fit$residuals
#
#Sample Size:                     10
#
#Test Statistic:                  W = 0.9613916
#
#Test Statistic Parameter:        n = 10
#
#P-value:                         0.801618
#
#Alternative Hypothesis:          True cdf does not equal the
#                                 Normal Distribution.


windows()
plot(TCE.fit, which = 1, 
	caption = "", 
	main = "Figure 21-5. Scatterplot of TCE Residuals vs. Regression Line Estimates")


# NOTE: typing the command plot(TCE.fit) will automatically produce
#       4 diagnostic plots involving the residuals



# Steps 3-4:  Compute a 95% upper confidence limit for TCE at months 30 and 26
#-----------------------------------------------------------------------------


# Confidence limits at Month = 30
#--------------------------------
predict(TCE.fit, newdata = data.frame(Month = 30), 
	interval = "confidence", level = 0.9)
#       fit      lwr      upr
#1 6.861126 2.380887 11.34136

# So the upper confidence limit associated with the 90% two-sided confidence limit
# (which is equivalent to the 95% upper confidence limit)
# is equal to 11.34 ppb, which is less than 20 ppb.


# NOTE: The Unified Guidance reports a value of UCL = 13.14, and the 
#       Errata update this value to 12.87.


# Confidence limits at Month = 26
#--------------------------------
predict(TCE.fit, newdata = data.frame(Month = 26), 
	interval = "confidence", level = 0.9)
#       fit    lwr      upr
#1 13.27170 9.6425 16.90091


# So the upper confidence limit associated with the 90% two-sided confidence limit
# (which is equivalent to the 95% upper confidence limit)
# is equal to 16.90 ppb, which is less than 20 ppb.

# NOTE: The Unified Guidance reports a value of UCL = 18.32, and the 
#       Errata update this value to 18.14.


# Step 5 - Compute upper 95% confidence limit ignoring linear trend
#------------------------------------------------------------------
enorm(TCE, ci = T, ci.type = "upper")
#
#Results of Distribution Parameter Estimation
#--------------------------------------------
#
#Assumed Distribution:            Normal
#
#Estimated Parameter(s):          mean = 30.42000
#                                 sd   = 15.50776
#
#Estimation Method:               mvue
#
#Data:                            TCE
#
#Sample Size:                     10
#
#Confidence Interval for:         mean
#
#Confidence Interval Method:      Exact
#
#Confidence Interval Type:        upper
#
#Confidence Level:                95%
#
#Confidence Interval:             LCL =     -Inf
#                                 UCL = 39.40956

rm(TCE, Month, TCE.fit)


#####################################################################################

# Example 21-8, pp. 21-31 to 21-32
#---------------------------------

EPA.09.Ex.17.7.sodium.df
#   Year Sodium.ppm
#1  89.6         56
#2  90.1         53
#3  90.8         51
#4  91.1         55
#5  92.1         52
#6  93.1         60
#7  94.1         62
#8  95.6         59
#9  96.1         61
#10 96.3         63


Na.kendall <- kendallTrendTest(Sodium.ppm ~ Year, 
	data = EPA.09.Ex.17.7.sodium.df, ci = T)
Na.kendall
#
#Results of Hypothesis Test
#--------------------------
#
#Null Hypothesis:                 tau = 0
#
#Alternative Hypothesis:          True tau is not equal to 0
#
#Test Name:                       Kendall's Test for Trend
#                                 (with continuity correction)
#
#Estimated Parameter(s):          tau       =   0.5555556
#                                 slope     =   1.3333333
#                                 intercept = -65.9666667
#
#Estimation Method:               slope:      Theil/Sen Estimator
#                                 intercept:  Conover's Estimator
#
#Data:                            y = Sodium.ppm
#                                 x = Year      
#
#Data Source:                     EPA.09.Ex.17.7.sodium.df
#
#Sample Size:                     10
#
#Test Statistic:                  z = 2.146625
#
#P-value:                         0.03182313
#
#Confidence Interval for:         slope
#
#Confidence Interval Method:      Gilbert's Modification
#                                 of Theil/Sen Method
#
#Confidence Interval Type:        two-sided
#
#Confidence Level:                95%
#
#Confidence Interval:             LCL = 0.3992083
#                                 UCL = 2.3333333

# NOTE:  Theil/Sen estimate of slope is 1.33 with a 
#        2-sided 95% CI of [0.40, 2.33]


# To reproduce Figure 21-6, do the following:

library(boot)

statistic.fcn <- function(data, indices, x.to.predict) {
	df <- data[indices, ]
	kendall.list <- kendallTrendTest(df$Sodium.ppm, df$Year)
	slope     <- kendall.list$estimate['slope']	
	intercept <- kendall.list$estimate['intercept']
	intercept + slope * x.to.predict	
}

n.points <- 101

Year <- EPA.09.Ex.17.7.sodium.df$Year
Na   <- EPA.09.Ex.17.7.sodium.df$Sodium.ppm

x.to.predict <- seq(min(Year), max(Year), length = n.points)

set.seed(27)
boot.out <- boot(data = EPA.09.Ex.17.7.sodium.df, 
	statistic = statistic.fcn, R = 500, 
	x.to.predict = x.to.predict)

LCL <- numeric(n.points)
UCL <- numeric(n.points)

for(i in 1:n.points) {
	boot.ci.obj <- boot.ci(boot.out, conf = 0.95, type = "perc", index = i)
	LCL[i] <- boot.ci.obj$percent[4]
	UCL[i] <- boot.ci.obj$percent[5]
}

windows()
plot(Year, Na, ylim = c(45, 70), pch = 16, 
	xlab = "Sample Year", ylab = "Sodium (ppm)",
	main = "Figure 21-6. 95% Theil-Sen Confidence Band on Sodium Measurements")
abline(a = Na.kendall$estimate['intercept'], b = Na.kendall$estimate['slope'])

lines(x.to.predict, LCL, lty = 3)
lines(x.to.predict, UCL, lty = 3)

legend(90, 70, c("Thiel-Sen Trend", "95% Conf Band"), lty = 1:2)


rm(Year, Na, Na.kendall, statistic.fcn, n.points, x.to.predict, 
	boot.out, LCL, UCL, i, boot.ci.obj)
