#####
# File:     Script - EPA09, Chapter 17 Examples.R
#
# Purpose:  Reproduce Examples in Chapter 17 of the EPA Guidance Document
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
# Updated:  October 03, 2012
#####

#NOTE:  Unless otherwise noted, differences between what is shown in the
#       EPA Guidance Document and results shown here are due to the 
#       EPA Guidance Document rounding intermediate results prior to 
#       the final result.

library(EnvStats)

#######################################################################################

# Example 17-1, pp. 17-6 to 17-8
#-------------------------------

# Raw data and summary statistics
#--------------------------------

head(EPA.09.Ex.17.1.loglead.df)
#  Month   Well  Well.type LogLead
#1     1 Well.1 Background    4.06
#2     2 Well.1 Background    3.99
#3     3 Well.1 Background    3.40
#4     4 Well.1 Background    3.83
#5     1 Well.2 Background    3.83
#6     2 Well.2 Background    4.34


longToWide(EPA.09.Ex.17.1.loglead.df, "LogLead", "Month", "Well", 
	paste.row.name = TRUE)
#        Well.1 Well.2 Well.3 Well.4 Well.5 Well.6
#Month.1   4.06   3.83   4.61   3.53   4.11   4.42
#Month.2   3.99   4.34   5.14   4.54   4.29   5.21
#Month.3   3.40   3.47   3.67   4.26   5.50   5.29
#Month.4   3.83   4.22   3.97   4.42   5.31   5.08


summaryStats(LogLead ~ Well, data = EPA.09.Ex.17.1.loglead.df,
	digits = 2, combine.groups = TRUE, 
	stats.in.rows = TRUE)
#       Well.1 Well.2 Well.3 Well.4 Well.5 Well.6 Combined
#N       4      4      4      4      4      4     24      
#Mean    3.82   3.96   4.35   4.19   4.8    5      4.35   
#SD      0.3    0.4    0.66   0.45   0.7    0.4    0.62   
#Median  3.91   4.03   4.29   4.34   4.8    5.14   4.28   
#Min     3.4    3.47   3.67   3.53   4.11   4.42   3.4    
#Max     4.06   4.34   5.14   4.54   5.5    5.29   5.5

summaryStats(LogLead ~ Well.type, data = EPA.09.Ex.17.1.loglead.df,
	digits = 2, combine.groups = TRUE, 
	stats.in.rows = TRUE)
#       Background Compliance Combined
#N       8         16         24      
#Mean    3.89       4.58       4.35   
#SD      0.33       0.61       0.62   
#Median  3.91       4.48       4.28   
#Min     3.4        3.53       3.4    
#Max     4.34       5.5        5.5



# Create new variable that combines background wells together but
# keeps compliance well separate
#----------------------------------------------------------------

new.Well <- as.character(EPA.09.Ex.17.1.loglead.df$Well)
new.Well[EPA.09.Ex.17.1.loglead.df$Well.type == "Background"] <- "Background"
new.Well <- factor(new.Well, levels = c("Background", paste("Well", 3:6, sep = ".")))



# Perform the ANOVA
#------------------

# Make sure Treatment contrasts are being used
options(contrasts = c("contr.treatment", "contr.poly"))

LogLead.fit <- lm(LogLead ~ new.Well, data = EPA.09.Ex.17.1.loglead.df)
anova(LogLead.fit)
#Analysis of Variance Table
#
#Response: LogLead
#          Df Sum Sq Mean Sq F value  Pr(>F)  
#new.Well   4 4.2888  1.0722  4.3852 0.01114 *
#Residuals 19 4.6456  0.2445                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

summary(LogLead.fit)
#Call:
#lm(formula = LogLead ~ new.Well)
#
#Residuals:
#     Min       1Q   Median       3Q      Max 
#-0.69250 -0.44000  0.08875  0.29937  0.79250 
#
#Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
#(Intercept)      3.8925     0.1748  22.265 4.49e-15 ***
#new.WellWell.3   0.4550     0.3028   1.503  0.14937    
#new.WellWell.4   0.2950     0.3028   0.974  0.34218    
#new.WellWell.5   0.9100     0.3028   3.005  0.00728 ** 
#new.WellWell.6   1.1075     0.3028   3.658  0.00167 ** 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
#
#Residual standard error: 0.4945 on 19 degrees of freedom
#Multiple R-squared:  0.48,      Adjusted R-squared: 0.3706 
#F-statistic: 4.385 on 4 and 19 DF,  p-value: 0.01114 

# The p-value for Well 3 is 0.15, for Well 4 is 0.34, so neither of 
# these is significant at the 0.05/4 = 0.0125 level.
# The p-value for Well 5 is 0.007 and the p-value for Well 6 is 0.002,
# both of which are less that 0.0125, so conclude that the means for 
# Wells 5 and 6 are significantly greater than the background mean.


rm(new.Well, LogLead.fit)


#######################################################################################

# Example 17-2, pp. 17-13 to 17-14
#---------------------------------

# Raw data and Ranks
#-------------------

head(EPA.09.Ex.17.2.toluene.df)
#  Month   Well  Well.type Toluene.ppb.orig Toluene.ppb Censored
#1     1 Well.1 Background               <5         5.0     TRUE
#2     2 Well.1 Background              7.5         7.5    FALSE
#3     3 Well.1 Background               <5         5.0     TRUE
#4     4 Well.1 Background               <5         5.0     TRUE
#5     5 Well.1 Background              6.4         6.4    FALSE
#6     1 Well.2 Background               <5         5.0     TRUE


longToWide(EPA.09.Ex.17.2.toluene.df, 
	"Toluene.ppb.orig", "Month", "Well", 
	paste.row.name = TRUE)
#        Well.1 Well.2 Well.3 Well.4 Well.5
#Month.1     <5     <5     <5     <5     <5
#Month.2    7.5     <5   12.5   13.7   20.1
#Month.3     <5     <5      8   15.3     35
#Month.4     <5     <5     <5   20.2   28.2
#Month.5    6.4     <5   11.2   25.1     19


longToWide(EPA.09.Ex.17.2.toluene.df, 
	"Toluene.ppb", "Month", "Well", 
	paste.row.name = TRUE)
#        Well.1 Well.2 Well.3 Well.4 Well.5
#Month.1    5.0      5    5.0    5.0    5.0
#Month.2    7.5      5   12.5   13.7   20.1
#Month.3    5.0      5    8.0   15.3   35.0
#Month.4    5.0      5    5.0   20.2   28.2
#Month.5    6.4      5   11.2   25.1   19.0



# Create new variable that combines background wells together but
# keep compliance well separate
#----------------------------------------------------------------

new.Well <- as.character(EPA.09.Ex.17.2.toluene.df$Well)
new.Well[EPA.09.Ex.17.2.toluene.df$Well.type == "Background"] <- "Background"
new.Well <- factor(new.Well, levels = c("Background", paste("Well", 3:5, sep = ".")))

summaryStats(rank(Toluene.ppb) ~ new.Well, 
	data = EPA.09.Ex.17.2.toluene.df, 
	stats.in.rows = TRUE, combine.groups = TRUE)[c("N", "Mean"), ]
#     Background Well.3 Well.4 Well.5 Combined
#N          10.0    5.0    5.0    5.0       25
#Mean        7.9   12.2   17.7   19.3       13



# Perform Kruskal-Wallis test
#----------------------------

kruskal.test(Toluene.ppb ~ new.Well, data = EPA.09.Ex.17.2.toluene.df)
#
#Results of Hypothesis Test
#--------------------------
#
#Alternative Hypothesis:          
#
#Test Name:                       Kruskal-Wallis rank sum test
#
#Data:                            Toluene by new.Well
#
#Test Statistic:                  Kruskal-Wallis chi-squared = 11.86932
#
#P-value:                         0.007844498



#Now do post-hoc comparisons
#---------------------------

# NOTE: Equation 17.14 on page 17-12 does not account for ties.
#       Here, we will perform parametric ANOVA on the ranks.

# Make sure Treatment contrasts are being used
options(contrasts = c("contr.treatment", "contr.poly"))

Toluene.fit <- lm(rank(Toluene.ppb) ~ new.Well, 
	data = EPA.09.Ex.17.2.toluene.df)

anova(Toluene.fit)
#Analysis of Variance Table
#
#Response: rank(Toluene)
#          Df Sum Sq Mean Sq F value   Pr(>F)   
#new.Well   3  572.2 190.733  6.8492 0.002148 **
#Residuals 21  584.8  27.848                    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


summary(Toluene.fit)
#
#Call:
#lm(formula = rank(Toluene) ~ new.Well, data = EPA.09.Ex.17.2.toluene.df)
#
#Residuals:
#   Min     1Q Median     3Q    Max 
# -12.8   -1.4    0.3    4.3    6.1 
#
#Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
#(Intercept)       7.900      1.669   4.734 0.000113 ***
#new.WellWell.3    4.300      2.890   1.488 0.151694    
#new.WellWell.4    9.800      2.890   3.391 0.002759 ** 
#new.WellWell.5   11.400      2.890   3.944 0.000742 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
#
#Residual standard error: 5.277 on 21 degrees of freedom
#Multiple R-squared: 0.4946,     Adjusted R-squared: 0.4223 
#F-statistic: 6.849 on 3 and 21 DF,  p-value: 0.002148 

# The p-value for Well 3 is 0.15, so it's not significant at the 0.05/3 = 0.0167 level.
# The p-value for Well 4 is 0.003 and the p-value for Well 5 is 0.0007,
# both of which are less that 0.0167, so conclude that the means for 
# Wells 4 and 5 are significantly greater than the background mean.


rm(new.Well, Toluene.fit)


#####################################################################################################

# Example 17-3, pp. 17-17 to 17-18
#---------------------------------

# Raw data and summary statistics
#--------------------------------

head(EPA.09.Ex.17.3.chrysene.df)
#  Month   Well  Well.type Chrysene.ppb
#1     1 Well.1 Background         19.7
#2     2 Well.1 Background         39.2
#3     3 Well.1 Background          7.8
#4     4 Well.1 Background         12.8
#5     1 Well.2 Background         10.2
#6     2 Well.2 Background          7.2

longToWide(EPA.09.Ex.17.3.chrysene.df, "Chrysene.ppb", "Month", "Well", 
	paste.row.name = TRUE)
#        Well.1 Well.2 Well.3 Well.4 Well.5
#Month.1   19.7   10.2   68.0   26.8   47.0
#Month.2   39.2    7.2   48.9   17.7   30.5
#Month.3    7.8   16.1   30.1   31.9   15.0
#Month.4   12.8    5.7   38.1   22.2   23.4


summaryStats(Chrysene.ppb ~ Well, data = EPA.09.Ex.17.3.chrysene.df,
	digits = 2, combine.groups = TRUE, stats.in.rows = TRUE)
#       Well.1 Well.2 Well.3 Well.4 Well.5 Combined
#N       4      4      4      4      4     20      
#Mean   19.88   9.8   46.27  24.65  28.98  25.91   
#SD     13.78   4.6   16.4    6.1   13.58  16.21   
#Median 16.25   8.7   43.5   24.5   26.95  22.8    
#Min     7.8    5.7   30.1   17.7   15      5.7    
#Max    39.2   16.1   68     31.9   47     68


# Shapiro-Wilk Test on Background observations
#---------------------------------------------

gofTest(Chrysene.ppb ~ 1, data = EPA.09.Ex.17.3.chrysene.df,
	subset = Well.type == "Background")
#
#Results of Goodness-of-Fit Test
#-------------------------------
#
#Test Method:                     Shapiro-Wilk GOF
#
#Hypothesized Distribution:       Normal
#
#Estimated Parameter(s):          mean = 14.83750
#                                 sd   = 10.92689
#
#Estimation Method:               mvue
#
#Data:                            Chrysene.ppb
#
#Subset With:                     Well.type == "Background"
#
#Data Source:                     EPA.09.Ex.17.3.chrysene.df
#
#Sample Size:                     8
#
#Test Statistic:                  W = 0.7978882
#
#Test Statistic Parameter:        n = 8
#
#P-value:                         0.02717253
#
#Alternative Hypothesis:          True cdf does not equal the
#                                 Normal Distribution.

 
# Shapiro-Wilk Test on Background observations, 
# assuming lognormal distribution
#----------------------------------------------

gofTest(Chrysene.ppb ~ 1, data = EPA.09.Ex.17.3.chrysene.df,
	subset = Well.type == "Background", 
	dist = "lnorm")
#
#Results of Goodness-of-Fit Test
#-------------------------------
#
#Test Method:                     Shapiro-Wilk GOF
#
#Hypothesized Distribution:       Lognormal
#
#Estimated Parameter(s):          meanlog = 2.5085773
#                                 sdlog   = 0.6279479
#
#Estimation Method:               mvue
#
#Data:                            Chrysene.ppb
#
#Subset With:                     Well.type == "Background"
#
#Data Source:                     EPA.09.Ex.17.3.chrysene.df
#
#Sample Size:                     8
#
#Test Statistic:                  W = 0.9564116
#
#Test Statistic Parameter:        n = 8
#
#P-value:                         0.7753094
#
#Alternative Hypothesis:          True cdf does not equal the
#                                 Lognormal Distribution.


# Raw data and summary statistics for log-transformed data
#---------------------------------------------------------

data.frame(Month = levels(Month), sapply(split(round(log(Chrysene), 3), Well), I))
#  Month Well.1 Well.2 Well.3 Well.4 Well.5
#1     1  2.981  2.322  4.220  3.288  3.850
#2     2  3.669  1.974  3.890  2.874  3.418
#3     3  2.054  2.779  3.405  3.463  2.708
#4     4  2.549  1.740  3.640  3.100  3.153


summaryStats(log(Chrysene.ppb) ~ Well, data = EPA.09.Ex.17.3.chrysene.df,
	digits = 3, combine.groups = TRUE, stats.in.rows = TRUE)
#       Well.1 Well.2 Well.3 Well.4 Well.5 Combined
#N       4      4      4      4      4     20      
#Mean    2.813  2.204  3.789  3.181  3.282  3.054  
#SD      0.685  0.452  0.349  0.253  0.479  0.681  
#Median  2.765  2.148  3.765  3.194  3.285  3.126  
#Min     2.054  1.74   3.405  2.874  2.708  1.74   
#Max     3.669  2.779  4.22   3.463  3.85   4.22


summaryStats(log(Chrysene.ppb) ~ Well.type, data = EPA.09.Ex.17.3.chrysene.df,
	digits = 3, combine.groups = TRUE, stats.in.rows = TRUE)
#       Background Compliance Combined
#N       8         12         20      
#Mean    2.509      3.417      3.054  
#SD      0.628      0.436      0.681  
#Median  2.436      3.411      3.126  
#Min     1.74       2.708      1.74   
#Max     3.669      4.22       4.22 


# Compute upper Tolerance Interval with 95% coverage and 95% confidence
# based on log-transformed Chrysene data for Background Wells
#----------------------------------------------------------------------

tolIntLnorm(Chrysene[Well.type == "Background"], coverage = 0.95, 
    ti.type = "upper", conf.level = 0.95)
#
#Results of Distribution Parameter Estimation
#--------------------------------------------
#
#Assumed Distribution:            Lognormal
#
#Estimated Parameter(s):          meanlog = 2.5085773
#                                 sdlog   = 0.6279479
#
#Estimation Method:               mvue
#
#Data:                            Chrysene[Well.type == "Background"]
#
#Sample Size:                     8
#
#Tolerance Interval Coverage:     95%
#
#Coverage Type:                   content
#
#Tolerance Interval Method:       Exact
#
#Tolerance Interval Type:         upper
#
#Confidence Level:                95%
#
#Tolerance Interval:              LTL =  0.0000
#                                 UTL = 90.9247


rm(Month, Well, Well.type, Chrysene)


#####################################################################################################

# Example 17-4, p. 17-21
#-----------------------

names(EPA.09.Ex.17.4.copper.df)
#[1] "Month"           "Well"            "Well.type"       "Copper.ppb.orig"
#[5] "Copper.ppb"      "Censored"   

Month     <- EPA.09.Ex.17.4.copper.df$Month
Well      <- EPA.09.Ex.17.4.copper.df$Well
Well.type <- EPA.09.Ex.17.4.copper.df$Well.type
Cu.orig   <- EPA.09.Ex.17.4.copper.df$Copper.ppb.orig
Cu        <- EPA.09.Ex.17.4.copper.df$Copper.ppb

data.frame(Month = levels(Month), sapply(split(Cu.orig, Well), I))
#  Month Well.1 Well.2 Well.3 Well.4 Well.5
#1     1     <5    9.2     <5              
#2     2     <5     <5    5.4              
#3     3    7.5     <5    6.7              
#4     4     <5    6.1     <5              
#5     5     <5      8     <5    6.2     <5
#6     6     <5    5.9     <5     <5     <5
#7     7    6.4     <5     <5    7.8    5.6
#8     8      6     <5     <5   10.4     <5


tolIntNpar(Cu[Well.type == "Background"], conf.level = 0.95, 
	ti.type = "upper", lb = 0) 
#
#Results of Distribution Parameter Estimation
#--------------------------------------------
#
#Assumed Distribution:            None
#
#Data:                            Cu[Well.type == "Background"]
#
#Sample Size:                     24
#
#Tolerance Interval Coverage:     88%
#
#Coverage Type:                   content
#
#Tolerance Interval Method:       Exact
#
#Tolerance Interval Type:         upper
#
#Confidence Level:                95%
#
#Tolerance Limit Rank(s):         24 
#
#Tolerance Interval:              LTL = 0.0
#                                 UTL = 9.2


rm(Month, Well, Well.type, Cu.orig, Cu)

#######################################################################################################

# Example 17-5, p. 17-26 to 17-30
#--------------------------------

names(EPA.09.Ex.17.5.chloride.df)
#[1] "Date"         "Chloride.ppm" "Elapsed.Days" "Residuals"   

Date   <- EPA.09.Ex.17.5.chloride.df$Date
Cl     <- EPA.09.Ex.17.5.chloride.df$Chloride.ppm


# Figure 17-2
#------------
windows()
plot(Date, Cl, pch = 16, 
	xlab = "Sampling Date", ylab = "Chloride (ppm)",
	main = paste("Figure 17-2. Time Series Plot of Chloride (ppm)",
	"Overlaid With Linear Regression", sep = "\n"))
Cl.fit <- lm(Chloride.ppm ~ Date, data = EPA.09.Ex.17.5.chloride.df)
abline(Cl.fit)


# Diagnostic Plots:  Figures 17-3 to 17-5
#----------------------------------------

windows()
qqPlot(Cl.fit$residuals, ylab = "Chloride Residuals (ppm)",
	main = "Figure 17-3. Probability Plot of \nChloride Regression Residuals")

windows()
plot(Date, Cl.fit$residuals, pch = 16,
	xlab = "Sampling Date", ylab = "Residual Chloride (ppm)",
	main = "Figure 17-4. Scatter Plot of \nChloride Residuals vs. Sampling Date")

windows()
plot(Cl.fit$fitted, Cl.fit$residuals, pch = 16,
	xlab = "Predicted Regression Line Fit (ppm)", 
	ylab = "Residual Chloride (ppm)",
	main = "Figure 17-5. Scatter Plot of \nChloride Residuals vs. Predicted Regression Fits")

# Alternatively, to automatically plot a subset of available diagnostic plots:
# ----------------------------------------------------------------------------
for(i in c(1:3, 5)){
	windows()
	plot(Cl.fit, which = i)
}


# Estimated regression coefficients and Associated p-values
#----------------------------------------------------------

summary(Cl.fit)$coef
#                 Estimate   Std. Error   t value     Pr(>|t|)
#(Intercept) -24.667125117 4.2285856469 -5.833422 1.994737e-05
#Date          0.003095753 0.0003345330  9.253952 4.764303e-08


rm(Date, Cl, Cl.fit, i)

############################################################################

# Example 17-6, p. 17-33 to 17-34
#--------------------------------

names(EPA.09.Ex.17.6.sulfate.df)
#[1] "Sample.No"     "Year"          "Month"         "Sampling.Date"
#[5] "Date"          "Sulfate.ppm" 

Sampling.Date <- EPA.09.Ex.17.6.sulfate.df$Sampling.Date
Sulfate <- EPA.09.Ex.17.6.sulfate.df$Sulfate.ppm

# Figure 17-6
#------------
windows()
plot(Sampling.Date, Sulfate, pch = 15, ylim = c(400, 900),
	xlab = "Sampling Date", ylab = "Sulfate Conc (ppm)",
	main = "Figure 17-6. Time Series Plot of \nSulfate Concentrations (ppm)")
Sulfate.fit <- lm(Sulfate.ppm ~ Sampling.Date, data = EPA.09.Ex.17.6.sulfate.df)
abline(Sulfate.fit, lty = 2)

# NOTE: Strictly speaking, Sampling.Date does not correctly represent time, 
#       since it has the month number in the 10'th place and there are 12
#       months, not 10.  A better plot is produced with the following code:

Date <- EPA.09.Ex.17.6.sulfate.df$Date

windows()
plot(Date, Sulfate, pch = 15, ylim = c(400, 900),
	xlab = "Sampling Date", ylab = "Sulfate Conc (ppm)",
	main = "Figure 17-6. Time Series Plot of \nSulfate Concentrations (ppm)")
Sulfate.fit <- lm(Sulfate.ppm ~ Date, data = EPA.09.Ex.17.6.sulfate.df)
abline(Sulfate.fit, lty = 2)



# Mann-Kendall Test for Trend using Sampling.Date
#------------------------------------------------

kendallTrendTest(y = Sulfate, x = Sampling.Date)
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
#Estimated Parameter(s):          tau       =     0.7667984
#                                 slope     =    26.6666667
#                                 intercept = -1909.3333333
#
#Estimation Method:               slope:      Theil/Sen Estimator
#                                 intercept:  Conover's Estimator
#
#Data:                            y = Sulfate      
#                                 x = Sampling.Date
#
#Sample Size:                     23
#
#Test Statistic:                  z = 5.107322
#
#P-value:                         3.267574e-07
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
#Confidence Interval:             LCL = 20.00000
#                                 UCL = 35.71182

rm(Sampling.Date, Sulfate, Sulfate.fit, Date)


#####################################################################################

# Example 17-7, p. 17-36 to 17-38
#--------------------------------

names(EPA.09.Ex.17.7.sodium.df)
#[1] "Year"       "Sodium.ppm"

Year <- EPA.09.Ex.17.7.sodium.df$Year
Na   <- EPA.09.Ex.17.7.sodium.df$Sodium.ppm

Kendall.list <- kendallTrendTest(Na, Year)
Kendall.list
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
#Data:                            y = Na  
#                                 x = Year
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

names(Kendall.list)
# [1] "statistic"         "parameters"        "p.value"          
# [4] "estimate"          "null.value"        "alternative"      
# [7] "method"            "estimation.method" "sample.size"      
#[10] "data.name"         "bad.obs"           "S"                
#[13] "var.S"             "slopes"            "interval"  

Kendall.list$estimate
#        tau       slope   intercept 
#  0.5555556   1.3333333 -65.9666667       


windows()
plot(Year, Na, ylim = c(50, 65), pch = 16, 
	xlab = "Year (2-Digit Format)", ylab = "Sodium (ppm)",
	main = "Figure 17-7. Time Series Plot of \nSodium Concentrations (ppm)")
abline(a = Kendall.list$estimate["intercept"], b = Kendall.list$estimate["slope"])
Na.fit <- lm(Sodium.ppm ~ Year, data = EPA.09.Ex.17.7.sodium.df)
abline(Na.fit, lty = 2)
legend("topleft", c("Thiel-Sen", "Linear Regression"), lty = c(1, 2))


rm(Year, Na, Kendall.list, Na.fit)

######################################################################################
