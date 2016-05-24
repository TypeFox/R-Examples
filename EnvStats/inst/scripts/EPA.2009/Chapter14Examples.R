#####
# File:     Script - EPA09, Chapter 14 Examples.R
#
# Purpose:  Reproduce Examples in Chapter 14 of the EPA Guidance Document
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
# Updated:  September 9, 2012
#####

#NOTE:  Unless otherwise noted, differences between what is shown in the
#       EPA Guidance Document and results shown here are due to the 
#       EPA Guidance Document rounding intermediate results prior to 
#       the final result.

#######################################################################################


# Example 14-1, pp. 14-5 to 14-6
#-------------------------------

head(EPA.09.Ex.14.1.manganese.df)
#  Quarter Well Manganese.ppm
#1       1 BW-1         28.14
#2       2 BW-1         29.33
#3       3 BW-1         30.45
#4       4 BW-1         32.42
#5       5 BW-1         34.37
#6       6 BW-1         33.25

Mn.df <- longToWide(EPA.09.Ex.14.1.manganese.df, 
	"Manganese.ppm", "Quarter", "Well",
	paste.row.name = TRUE)
Mn.df
#           BW-1  BW-2  BW-3  BW-4
#Quarter.1 28.14 31.41 27.15 30.46
#Quarter.2 29.33 30.27 30.24 30.60
#Quarter.3 30.45 32.57 29.14 30.96
#Quarter.4 32.42 32.77 30.59 30.70
#Quarter.5 34.37 33.03 34.88 32.71
#Quarter.6 33.25 32.18 30.53 31.76
#Quarter.7 31.02 28.85 30.33 31.85
#Quarter.8 28.50 32.88 30.42 29.58

n.Wells <- ncol(Mn.df)
windows()
matplot(y = Mn.df, type = "l", lty = 1:n.Wells, col = 1:n.Wells, lwd = 2, 
	xlim = c(0, 10), ylim = c(25, 40), 
	xlab = "Quarter", ylab = "Manganese (ppm)",
	main = "Figure 14-2. Manganese Parallel Time Series Plot")
legend(0.5, 40, names(Mn.df), lty = 1:n.Wells, 
	lwd = rep(2, n.Wells), col = 1:n.Wells)

rm(Mn.df, n.Wells)


#######################################################################################


# Example 14-2, pp. 14-9 to 14-12
#--------------------------------

# NOTE: these data are balanced (i.e., equal number of quarters per well).

# Step 1
#-------

Mn.df <- longToWide(EPA.09.Ex.14.1.manganese.df, 
	"Manganese.ppm", "Quarter", "Well",
	paste.row.name = TRUE)
Mn.df
#           BW-1  BW-2  BW-3  BW-4
#Quarter.1 28.14 31.41 27.15 30.46
#Quarter.2 29.33 30.27 30.24 30.60
#Quarter.3 30.45 32.57 29.14 30.96
#Quarter.4 32.42 32.77 30.59 30.70
#Quarter.5 34.37 33.03 34.88 32.71
#Quarter.6 33.25 32.18 30.53 31.76
#Quarter.7 31.02 28.85 30.33 31.85
#Quarter.8 28.50 32.88 30.42 29.58

Mn.w.Qtr.Means.df <- data.frame(Event.Mean = apply(Mn.df, 1, mean), Mn.df)
Mn.w.Qtr.Means.df
#          Event.Mean  BW.1  BW.2  BW.3  BW.4
#Quarter.1    29.2900 28.14 31.41 27.15 30.46
#Quarter.2    30.1100 29.33 30.27 30.24 30.60
#Quarter.3    30.7800 30.45 32.57 29.14 30.96
#Quarter.4    31.6200 32.42 32.77 30.59 30.70
#Quarter.5    33.7475 34.37 33.03 34.88 32.71
#Quarter.6    31.9300 33.25 32.18 30.53 31.76
#Quarter.7    30.5125 31.02 28.85 30.33 31.85
#Quarter.8    30.3450 28.50 32.88 30.42 29.58

Grand.Mean <- mean(unlist(Mn.df))
Grand.Mean
#[1] 31.04188


# Step 2
#-------

fit <- aov(Manganese.ppm ~ factor(Quarter), 
	data = EPA.09.Ex.14.1.manganese.df)
Residuals <- fit$residuals

new.df <- data.frame(EPA.09.Ex.14.1.manganese.df, Residuals)

longToWide(new.df, "Residuals", "Quarter", "Well", 
	paste.row.name = TRUE)
#             BW-1    BW-2    BW-3    BW-4
#Quarter.1 -1.1500  2.1200 -2.1400  1.1700
#Quarter.2 -0.7800  0.1600  0.1300  0.4900
#Quarter.3 -0.3300  1.7900 -1.6400  0.1800
#Quarter.4  0.8000  1.1500 -1.0300 -0.9200
#Quarter.5  0.6225 -0.7175  1.1325 -1.0375
#Quarter.6  1.3200  0.2500 -1.4000 -0.1700
#Quarter.7  0.5075 -1.6625 -0.1825  1.3375
#Quarter.8 -1.8450  2.5350  0.0750 -0.7650


# Step 3
#-------

windows()
qqPlot(Residuals, add.line = T,
	ylab = "Event Mean Residuals (ppm)",
	main = "Figure 14-3. Probability Plot of Manganese Sampling Event Residuals",
	cex.main = 1.1)

gofTest(Residuals, test = "ppcc")

#Results of Goodness-of-Fit Test
#-------------------------------
#
#Test Method:                     PPCC GOF
#
#Hypothesized Distribution:       Normal
#
#Estimated Parameter(s):          mean = -3.989864e-17
#                                 sd   =  1.203236e+00
#
#Estimation Method:               mvue
#
#Data:                            Residuals
#
#Sample Size:                     32
#
#Test Statistic:                  r = 0.9935174
#
#Test Statistic Parameter:        n = 32
#
#P-value:                         0.9125129
#
#Alternative Hypothesis:          True cdf does not equal the
#                                 Normal Distribution.


# Step 4
#-------
varGroupTest(Residuals ~ Quarter, data = new.df)

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
#Estimated Parameter(s):          Quarter.1 = 3.921800
#                                 Quarter.2 = 0.297000
#                                 Quarter.3 = 2.011667
#                                 Quarter.4 = 1.289933
#                                 Quarter.5 = 1.087092
#                                 Quarter.6 = 1.264600
#                                 Quarter.7 = 1.614558
#                                 Quarter.8 = 3.473700
#
#Data:                            Residuals
#
#Grouping Variable:               Quarter
#
#Data Source:                     new.df
#
#Sample Sizes:                    Quarter.1 = 4
#                                 Quarter.2 = 4
#                                 Quarter.3 = 4
#                                 Quarter.4 = 4
#                                 Quarter.5 = 4
#                                 Quarter.6 = 4
#                                 Quarter.7 = 4
#                                 Quarter.8 = 4
#
#Test Statistic:                  F = 1.299012
#
#Test Statistic Parameters:       num df   =  7
#                                 denom df = 24
#
#P-value:                         0.2930252


# Steps 5-7
#----------

anova.fit <- anova(fit)
anova.fit
#Analysis of Variance Table
#
#Response: Mn
#                Df Sum Sq Mean Sq F value   Pr(>F)   
#factor(Quarter)  7 52.861  7.5516  4.0382 0.004685 **
#Residuals       24 44.881  1.8700                    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 


# Step 8
#-------

n.Wells <- ncol(Mn.df)
n.Qtrs  <- nrow(Mn.df)

Adj.SD <- sqrt( (1/n.Wells) * ( 
	anova.fit["factor(Quarter)", "Mean Sq"] + 
		(n.Wells - 1) * anova.fit["Residuals", "Mean Sq"] ) )
Adj.SD
#[1] 1.813955

#NOTE:  Equation [14.13] for the effective sample size is updated in the
#       EPA Errata Sheet.  All occurrences of T are replaced with TK.

F.T <- anova.fit["factor(Quarter)", "F value"]
Effective.n <- 1 + ( 
	( n.Qtrs * (n.Qtrs - 1) * (F.T + n.Wells - 1)^2 ) / 
	( n.Qtrs * F.T^2 + (n.Qtrs - 1) * (n.Wells - 1) )
	)

Effective.n
#[1] 19.31570


rm(Mn.df, Mn.w.Qtr.Means.df, Grand.Mean, 
	fit, Residuals, new.df, 
	anova.fit, n.Wells, n.Qtrs, Adj.SD, F.T, Effective.n)


#####################################################################################


# Example 14-3, pp. 14-14 to 14-16
#---------------------------------

# Step 1
#-------

head(EPA.09.Ex.14.3.alkalinity.df)
#        Date Alkalinity.mg.per.L
#1 1996-01-26                1400
#2 1996-02-20                1700
#3 1996-03-19                1900
#4 1996-04-22                1800
#5 1996-05-22                1300
#6 1996-06-24                2000

Alk  <- EPA.09.Ex.14.3.alkalinity.df$Alkalinity.mg.per.L
Date <- EPA.09.Ex.14.3.alkalinity.df$Date

windows()
plot(Date, Alk, type = "b", pch = 16,
	xlab = "Sampling Date", ylab = "Total Alkalinity (mg/L)",
	main = "Figure 14-4. Time Series Plot of Total Alkalinity (mg/L)")


# Steps 2-5
#----------
windows()
acf.list <- acf(Alk, 
	main = "Figure 14-5. Sample Autocorrelation Function for Total Alkalinity",
	cex.main = 1.1)
acf.list
#
#Autocorrelations of series ‘Alk’, by lag
#
#     0      1      2      3      4      5      6      7      8      9     10 
# 1.000  0.641  0.313  0.043 -0.230 -0.453 -0.496 -0.374 -0.174  0.040  0.171 
#    11     12     13     14     15     16     17 
# 0.303  0.416  0.368  0.228 -0.015 -0.175 -0.382

rm(Alk, Date, acf.list)


######################################################################################


# Example 14-4, p. 14-18
#-----------------------

serialCorrelationTest(EPA.09.Ex.14.4.arsenic.df$Arsenic.ppb)
#
#Results of Hypothesis Test
#--------------------------
#
#Null Hypothesis:                 rho = 0
#
#Alternative Hypothesis:          True rho is not equal to 0
#
#Test Name:                       Rank von Neumann Test for
#                                 Lag-1 Autocorrelation
#                                 (Beta Approximation)
#
#Estimated Parameter(s):          rho = 0.1084613
#
#Estimation Method:               Yule-Walker
#
#Data:                            EPA.09.Ex.14.4.arsenic.df$Arsenic.ppb
#
#Sample Size:                     16
#
#Test Statistic:                  RVN = 1.670588
#
#P-value:                         0.5051189
#
#Confidence Interval for:         rho
#
#Confidence Interval Method:      Normal Approximation
#
#Confidence Interval Type:        two-sided
#
#Confidence Level:                95%
#
#Confidence Interval:             LCL = -0.3786391
#                                 UCL =  0.5955617


####################################################################################


# Example 14-8, p. 14-32 to 14-33
#--------------------------------


# NOTE: The Adjusted Concentrations are already part of EPA.09.Ex.14.8.df.
#       In this example, we will recreate them.

head(EPA.09.Ex.14.8.df)
#     Month Year Unadj.Conc Adj.Conc
#1  January 1983       1.99     2.11
#2 February 1983       2.10     2.14
#3    March 1983       2.12     2.10
#4    April 1983       2.12     2.13
#5      May 1983       2.11     2.12
#6     June 1983       2.15     2.12

Unadj.df <- longToWide(EPA.09.Ex.14.8.df, "Unadj.Conc", "Month", "Year")
Unadj.df
#          1983 1984 1985
#January   1.99 2.01 2.15
#February  2.10 2.10 2.17
#March     2.12 2.17 2.27
#April     2.12 2.13 2.23
#May       2.11 2.13 2.24
#June      2.15 2.18 2.26
#July      2.19 2.25 2.31
#August    2.18 2.24 2.32
#September 2.16 2.22 2.28
#October   2.08 2.13 2.22
#November  2.05 2.08 2.19
#December  2.08 2.16 2.22

Monthly.Avg <- apply(Unadj.df, 1, mean)
Three.Yr.Mean <- mean(EPA.09.Ex.14.8.df$Unadj.Conc)
Adj.df <- round(Unadj.df - Monthly.Avg + Three.Yr.Mean, 2)

names(Unadj.df) <- paste("Unadj", names(Unadj.df), sep = ".")
names(Adj.df) <- paste("Adj", names(Adj.df), sep = ".")

EPA.df <- data.frame(Unadj.df, Monthly.Avg = round(Monthly.Avg, 2), Adj.df)
EPA.df
#          Unadj.1983 Unadj.1984 Unadj.1985 Monthly.Avg Adj.1983 Adj.1984 Adj.1985
#January         1.99       2.01       2.15        2.05     2.11     2.13     2.27
#February        2.10       2.10       2.17        2.12     2.14     2.14     2.21
#March           2.12       2.17       2.27        2.19     2.10     2.15     2.25
#April           2.12       2.13       2.23        2.16     2.13     2.14     2.24
#May             2.11       2.13       2.24        2.16     2.12     2.14     2.25
#June            2.15       2.18       2.26        2.20     2.12     2.15     2.23
#July            2.19       2.25       2.31        2.25     2.11     2.17     2.23
#August          2.18       2.24       2.32        2.25     2.10     2.16     2.24
#September       2.16       2.22       2.28        2.22     2.11     2.17     2.23
#October         2.08       2.13       2.22        2.14     2.10     2.15     2.24
#November        2.05       2.08       2.19        2.11     2.11     2.14     2.25
#December        2.08       2.16       2.22        2.15     2.09     2.17     2.23

round(Three.Yr.Mean, 2)
#[1] 2.17


# Now plot the data

Unadj.Conc <- EPA.09.Ex.14.8.df$Unadj.Conc
Adj.Conc   <- EPA.09.Ex.14.8.df$Adj.Conc
n          <- nrow(EPA.09.Ex.14.8.df)
Time <- paste(substr(EPA.09.Ex.14.8.df$Month, 1, 3), 
	substr(EPA.09.Ex.14.8.df$Year, 3, 4), sep = "-")

windows()
par(mar = c(7, 4, 3, 1) + 0.1, cex.lab = 1.25)
plot(1:n, Unadj.Conc, type = "n", xaxt = "n",
	xlab = "Time (Month)", 
	ylab = "ANALYTE CONCENTRATION (mg/L)", 
	main = "Figure 14-15. Seasonal Time Series Over a Three Year Period",
	cex.main = 1.1)
axis(1, at = 1:n, labels = rep("", n))
at <- (1:n)[substr(Time, 1, 3) %in% c("Jan", "May", "Sep")]
axis(1, at = at, labels = Time[at], cex.axis = 0.8)
points(1:n, Unadj.Conc, pch = 0, type = "o", lwd = 2)
points(1:n, Adj.Conc, pch = 3, type = "o", col = 8, lwd = 2)
abline(h = Three.Yr.Mean, lwd = 2)
legend("topleft", c("Unadjusted", "Adjusted", "3-Year Mean"), bty = "n", 
	pch = c(0, 3, -1), lty = c(1, 1, 1), lwd = 2, col = c(1, 8, 1),
	inset = c(0.05, 0.01))

rm(Unadj.df, Monthly.Avg, Three.Yr.Mean, Adj.df, EPA.df, 
	Unadj.Conc, Adj.Conc, n, Time, at)


#####################################################################################


# Example 14-9, p. 14-35 to 14-36
#--------------------------------

# Step 1
#-------
Mn.df <- longToWide(EPA.09.Ex.14.1.manganese.df, 
	"Manganese.ppm", "Quarter", "Well",
	paste.row.name = TRUE)
Mn.df
#           BW-1  BW-2  BW-3  BW-4
#Quarter.1 28.14 31.41 27.15 30.46
#Quarter.2 29.33 30.27 30.24 30.60
#Quarter.3 30.45 32.57 29.14 30.96
#Quarter.4 32.42 32.77 30.59 30.70
#Quarter.5 34.37 33.03 34.88 32.71
#Quarter.6 33.25 32.18 30.53 31.76
#Quarter.7 31.02 28.85 30.33 31.85
#Quarter.8 28.50 32.88 30.42 29.58

Mn.w.Qtr.Means.df <- data.frame(Event.Mean = apply(Mn.df, 1, mean), Mn.df)
Mn.w.Qtr.Means.df
#          Event.Mean  BW.1  BW.2  BW.3  BW.4
#Quarter.1    29.2900 28.14 31.41 27.15 30.46
#Quarter.2    30.1100 29.33 30.27 30.24 30.60
#Quarter.3    30.7800 30.45 32.57 29.14 30.96
#Quarter.4    31.6200 32.42 32.77 30.59 30.70
#Quarter.5    33.7475 34.37 33.03 34.88 32.71
#Quarter.6    31.9300 33.25 32.18 30.53 31.76
#Quarter.7    30.5125 31.02 28.85 30.33 31.85
#Quarter.8    30.3450 28.50 32.88 30.42 29.58

Grand.Mean <- mean(unlist(Mn.df))
Grand.Mean
#[1] 31.04188

fit <- aov(Manganese.ppm ~ factor(Quarter), 
	data = EPA.09.Ex.14.1.manganese.df)
Residuals <- fit$residuals

new.EPA.df <- data.frame(EPA.09.Ex.14.1.manganese.df, Residuals)

Step.1.df <- data.frame(Event.Mean = Mn.w.Qtr.Means.df$Event.Mean, 
	longToWide(new.EPA.df, "Residuals", "Quarter", "Well", 
	paste.row.name = TRUE), check.names = FALSE)
Step.1.df
#          Event.Mean    BW-1    BW-2    BW-3    BW-4
#Quarter.1    29.2900 -1.1500  2.1200 -2.1400  1.1700
#Quarter.2    30.1100 -0.7800  0.1600  0.1300  0.4900
#Quarter.3    30.7800 -0.3300  1.7900 -1.6400  0.1800
#Quarter.4    31.6200  0.8000  1.1500 -1.0300 -0.9200
#Quarter.5    33.7475  0.6225 -0.7175  1.1325 -1.0375
#Quarter.6    31.9300  1.3200  0.2500 -1.4000 -0.1700
#Quarter.7    30.5125  0.5075 -1.6625 -0.1825  1.3375
#Quarter.8    30.3450 -1.8450  2.5350  0.0750 -0.7650


# Step 2
#-------
Step.2.df <- Step.1.df
Step.2.df[, "Event.Mean"] <- round(Step.2.df[, "Event.Mean"], 3)
Step.2.df[, -1] <- round(Step.2.df[, -1] + Grand.Mean, 2)
Step.2.df
#          Event.Mean  BW.1  BW.2  BW.3  BW.4
#Quarter.1     29.290 29.89 33.16 28.90 32.21
#Quarter.2     30.110 30.26 31.20 31.17 31.53
#Quarter.3     30.780 30.71 32.83 29.40 31.22
#Quarter.4     31.620 31.84 32.19 30.01 30.12
#Quarter.5     33.748 31.66 30.32 32.17 30.00
#Quarter.6     31.930 32.36 31.29 29.64 30.87
#Quarter.7     30.512 31.55 29.38 30.86 32.38
#Quarter.8     30.345 29.20 33.58 31.12 30.28


# Step 3
#-------
n.Wells <- ncol(Step.2.df)
windows()
matplot(y = Step.2.df[, -1], type = "l", 
	lty = 1:n.Wells, col = 1:n.Wells, lwd = 2, 
	xlim = c(0, 10), ylim = c(28, 36), 
	xlab = "Quarter", ylab = "Adjusted Manganese (ppm)",
	main = "Figure 14-16. Parallel Time Series Plot of Adjusted Manganese Concentrations", 
	cex.main = 1.1)
legend("topleft", names(Step.2.df[, -1]), lty = 1:n.Wells, 
	lwd = rep(2, n.Wells), col = 1:n.Wells, bty = "n")

rm(Mn.df, Mn.w.Qtr.Means.df, Grand.Mean, fit, Residuals, 
	new.EPA.df, Step.1.df, Step.2.df, n.Wells)


#####################################################################################


# Example 14-10, p. 14-38 to 14-40
#---------------------------------

kendallSeasonalTrendTest(Unadj.Conc ~ Month + Year, 
	data = EPA.09.Ex.14.8.df)


#Results of Hypothesis Test
#--------------------------
#
#Null Hypothesis:                 All 12 values of tau = 0
#
#Alternative Hypothesis:          The seasonal taus are not all equal
#                                 (Chi-Square Heterogeneity Test)
#                                 At least one seasonal tau != 0
#                                 and all non-zero tau's have the
#                                 same sign (z Trend Test)
#
#Test Name:                       Seasonal Kendall Test for Trend
#                                 (with continuity correction)
#
#Estimated Parameter(s):          tau       =    0.9722222
#                                 slope     =    0.0600000
#                                 intercept = -131.7350000
#
#Estimation Method:               tau:        Weighted Average of
#                                             Seasonal Estimates
#                                 slope:      Hirsch et al.'s
#                                             Modification of
#                                             Thiel/Sen Estimator
#                                 intercept:  Median of
#                                             Seasonal Estimates
#
#Data:                            y      = Unadj.Conc
#                                 season = Month     
#                                 year   = Year      
#
#Data Source:                     EPA.09.Ex.14.8.df
#
#Sample Sizes:                    January   =  3
#                                 February  =  3
#                                 March     =  3
#                                 April     =  3
#                                 May       =  3
#                                 June      =  3
#                                 July      =  3
#                                 August    =  3
#                                 September =  3
#                                 October   =  3
#                                 November  =  3
#                                 December  =  3
#                                 Total     = 36
#
#Test Statistics:                 Chi-Square (Het) = 0.1071882
#                                 z (Trend)        = 5.1849514
#
#Test Statistic Parameter:        df = 11
#
#P-values:                        Chi-Square (Het) = 1.000000e+00
#                                 z (Trend)        = 2.160712e-07
#
#Confidence Interval for:         slope
#
#Confidence Interval Method:      Gilbert's Modification of
#                                 Theil/Sen Method
#
#Confidence Interval Type:        two-sided
#
#Confidence Level:                95%
#
#Confidence Interval:             LCL = 0.05786914
#                                 UCL = 0.07213086
