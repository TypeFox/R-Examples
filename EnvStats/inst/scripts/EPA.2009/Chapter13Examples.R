#####
# File:     Script - EPA09, Chapter 13 Examples.R
#
# Purpose:  Reproduce Examples in Chapter 13 of the EPA Guidance Document
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
# Updated:  September 3, 2012
#####

#NOTE:  Unless otherwise noted, differences between what is shown in the
#       EPA Guidance Document and results shown here are due to the 
#       EPA Guidance Document rounding intermediate results prior to 
#       the final result.

#######################################################################################

library(EnvStats)


# Example 13-1, pp. 13-3 to 13-5
#-------------------------------

#NOTE:  The EPA Errata Sheet shows updated values for some summary statistics.

# Stats for Iron, by Well

stats <- summaryStats(Iron.ppm ~ Well, data = EPA.09.Ex.13.1.iron.df, 
	combine.groups = FALSE, digits = 2, stats.in.rows = TRUE)
stats
#       Well.1 Well.2 Well.3 Well.4 Well.5 Well.6
#N        4      4      4      4      4      4   
#Mean    47.01  55.73  90.86  70.43 145.24 156.32
#SD      12.4   20.34  59.35  25.95  92.16  51.2 
#Median  50.05  57.05  76.73  76.95 137.66 171.93
#Min     29.96  32.14  39.25  34.12  60.95  83.1 
#Max     57.97  76.71 170.72  93.69 244.69 198.34


#####
# Figure 13-1

windows()
boxplot(Iron.ppm ~ Well, data = EPA.09.Ex.13.1.iron.df, 
	ylim = c(0, 250), ylab = "Iron (ppm)",
	main = "Figure 13-1. Side-by-Side Iron Boxplots")
points(1:ncol(stats), stats["Mean", ], pch = 5)

# NOTE: with so few observations per well, it is often useful to 
#       look at the raw data (optionally with confidence intervals)

windows()
stripChart(Iron.ppm ~ Well, data = EPA.09.Ex.13.1.iron.df, 
	show.ci = F, location.scale.text.cex = 0.7,
	ylim = c(0, 250), ylab = "Iron (ppm)")
mtext("Alternative Figure 13-1. Side-by-Side Iron Strip Plots", 
	cex = 1.25, font = 2, line = 2.5)

windows()
stripChart(Iron.ppm ~ Well, data = EPA.09.Ex.13.1.iron.df,
	location.scale.text.cex = 0.7,
	ylim = c(0, 300), ylab = "Iron (ppm)")
mtext(paste("Alternative Figure 13-1. Side-by-Side Iron Strip Plots",
	"with 95% Confidence Intervals", sep = "\n"), 
	cex = 1.25, font = 2, line = 1.75)


# Stats for log(Iron), by Well

stats <- summaryStats(log(Iron.ppm) ~ Well, data = EPA.09.Ex.13.1.iron.df, 
	combine.groups = F, digits = 2, stats.in.rows = TRUE)
stats
#       Well.1 Well.2 Well.3 Well.4 Well.5 Well.6
#N      4      4      4      4      4      4     
#Mean   3.82   3.97   4.35   4.19   4.8    5     
#SD     0.3    0.4    0.66   0.45   0.7    0.4   
#Median 3.91   4.02   4.29   4.34   4.8    5.14  
#Min    3.4    3.47   3.67   3.53   4.11   4.42  
#Max    4.06   4.34   5.14   4.54   5.5    5.29

#####
# Figure 13-2

windows()
boxplot(log(Iron.ppm) ~ Well, data = EPA.09.Ex.13.1.iron.df, 
	ylim = c(3, 6), ylab = "Log Iron log(ppm)",
	main = "Figure 13-2. Side-by-Side Log(Iron) Boxplots")
points(1:ncol(stats), stats["Mean",], pch = 5)


windows()
stripChart(log(Iron.ppm) ~ Well, data = EPA.09.Ex.13.1.iron.df, 
	show.ci = F, location.scale.text.cex = 0.7,
	ylim = c(3, 6), ylab = "Log Iron log(ppm)")
mtext("Alternative Figure 13-2. Side-by-Side Log(Iron) Strip Plots",
	cex = 1.15, font = 2, line = 2.5)

windows()
stripChart(log(Iron.ppm) ~ Well, data = EPA.09.Ex.13.1.iron.df,
	location.scale.text.cex = 0.7, 
	ylim = c(3, 6), ylab = "Log(Iron) log(ppm)")
mtext(paste("Alternative Figure 13-2. Side-by-Side Log(Iron) Strip Plots",
	"with 95% Confidence Intervals", sep = "\n"), line = 1.75, 
	cex = 1.15, font = 2)


rm(stats)


#######################################################################################


# Example 13-2, pp. 13-7 to 13-8
#-------------------------------

#NOTE: The EPA Errata Sheet gives updated values for this example.

stats <- summaryStats(log(Iron.ppm) ~ Well, data = EPA.09.Ex.13.1.iron.df,
	digits = 3, combine.groups = TRUE, stats.in.rows = TRUE)
stats
#       Well.1 Well.2 Well.3 Well.4 Well.5 Well.6 Combined
#N       4      4      4      4      4      4     24      
#Mean    3.82   3.965  4.347  4.187  4.803  5      4.354  
#SD      0.296  0.395  0.658  0.453  0.704  0.396  0.623  
#Median  3.91   4.025  4.29   4.34   4.8    5.145  4.275  
#Min     3.4    3.47   3.67   3.53   4.11   4.42   3.4    
#Max     4.06   4.34   5.14   4.54   5.5    5.29   5.5 

fit <- aov(log(Iron.ppm) ~ Well, data = EPA.09.Ex.13.1.iron.df)
anova(fit)
#Analysis of Variance Table
#
#Response: log(Iron.ppm)
#          Df Sum Sq Mean Sq F value  Pr(>F)  
#Well       5 4.3313 0.86626  3.3871 0.02486 *
#Residuals 18 4.6036 0.25575                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 


rm(fit, stats)


########################################################################################


# Example 13-3, pp. 13-11 to 13-10
#---------------------------------

#NOTE:  The EPA Errata Sheet gives updated values for this example.

stats <- summaryStats(log(Iron.ppm) ~ Well, 
	data = EPA.09.Ex.13.1.iron.df, 
	combine.groups = F, digits = 3, 
	stats.in.rows = TRUE)

# Unadjusted 99% Upper Prediction Limits
#---------------------------------------
UPLs <- sapply(
	with(EPA.09.Ex.13.1.iron.df, split(Iron.ppm, Well)), 
	function(x, pi.type, conf.level) {
		predIntLnorm(x, pi.type = "upper", 
			conf.level = 0.99)$interval$limits["UPL"]
	})
rbind(stats[c("Mean", "SD", "N"), ], "99% PL" = round(UPLs, 1))
#        Well.1  Well.2   Well.3  Well.4   Well.5   Well.6
#Mean     3.820   3.965    4.347   4.187    4.803    5.000
#SD       0.296   0.395    0.658   0.453    0.704    0.396
#N        4.000   4.000    4.000   4.000    4.000    4.000
#99% PL 205.000 392.200 2181.000 657.200 4340.000 1109.200



# Adjusted 99% Upper Prediction Limits
#-------------------------------------
fit <- aov(log(Iron.ppm) ~ Well, data = EPA.09.Ex.13.1.iron.df)
anova.fit <- anova(fit)
anova.fit
#Analysis of Variance Table
#
#Response: log(Iron.ppm)
#          Df Sum Sq Mean Sq F value  Pr(>F)  
#Well       5 4.3313 0.86626  3.3871 0.02486 *
#Residuals 18 4.6036 0.25575                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

RMSE <- anova.fit["Residuals", "Mean Sq"]
Df  <- anova.fit["Residuals", "Df"]
t.crit <- qt(0.99, df = Df)

UPLs <- sapply(
	with(EPA.09.Ex.13.1.iron.df, split(log(Iron.ppm), Well)), 
	function(x, t.crit, RMSE) {
		exp(mean(x) + t.crit * sqrt(RMSE * (1 + 1/length(x))))
	}, 
	t.crit = t.crit, RMSE = RMSE)
rbind(stats[c("Mean", "SD", "N"), ], "99% PL" = round(UPLs, 1))
#        Well.1  Well.2  Well.3  Well.4  Well.5  Well.6
#Mean     3.820   3.965   4.347   4.187   4.803   5.000
#SD       0.296   0.395   0.658   0.453   0.704   0.396
#N        4.000   4.000   4.000   4.000   4.000   4.000
#99% PL 193.100 223.200 327.200 278.800 515.800 628.400

rm(stats, UPLs, fit, anova.fit, RMSE, Df, t.crit)


