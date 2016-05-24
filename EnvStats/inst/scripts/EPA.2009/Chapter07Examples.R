#####
# File:     Script - EPA09, Chapter 07 Examples.R
#
# Purpose:  Reproduce Examples in Chapter 7 of the EPA Guidance Document
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
# Updated:  August 20, 2012
#####

#NOTE:  Unless otherwise noted, differences between what is shown in the
#       EPA Guidance Document and results shown here are due to the 
#       EPA Guidance Document rounding intermediate results prior to 
#       the final result.

#######################################################################################

library(EnvStats)

# Example 7-1, pp. 7-26 to 7-28
#------------------------------

stats <- summaryStats(Arsenic.ug.per.l ~ Data.Source, 
	data = EPA.09.Ex.7.1.arsenic.df,
	combine.groups = F, digits = 2, stats.in.rows = TRUE)
stats
#       Historical Case.1 Case.2
#N       8          4      4    
#Mean   37         69.58  82.92 
#SD     18.16      11.15  11.24 
#Median 36.4       69.15  81.7  
#Min    10.8       58.7   73.3  
#Max    74.1       81.3   95 


#####
# 1-of-2 future values option

kappa <- predIntNormSimultaneousK(n = 8, k = 1, m = 2, r = 6, 
	rule = "k.of.m", pi.type = "upper", 
	conf.level = 0.95)
kappa <- round(kappa, 2)
kappa
#[1] 1.8
# NOTE:  EPA Guidance Document states kappa = 1.83


Arsenic <- EPA.09.Ex.7.1.arsenic.df$Arsenic.ug.per.l
Source <- EPA.09.Ex.7.1.arsenic.df$Data.Source

UPL <- predIntNormSimultaneous(x = Arsenic[Source == "Historical"], 
	k = 1, m = 2, r = 6, 
	rule="k.of.m", pi.type="upper", conf.level = 0.95)$interval$limits["UPL"]
UPL <- round(UPL, 2)
UPL
# UPL 
#69.7
# NOTE:  EPA Guidance Document states UPL = 70.2 ug/l


#####
# Upper Tolerance Limit
UTL <- tolIntNorm(x = Arsenic[Source == "Historical"], 
	coverage = 0.95, cov.type = "content", 
	ti.type = "upper", conf.level = 0.95)$interval$limits["UTL"]

UTL <- round(UTL, 1)
UTL
# UTL 
#94.9 



#####
# 90% and 95% Confidence Limits for the Mean for Case 1 and Case 2

# Case 1
ci.lower.limit.90.Pct <- enorm(Arsenic[Source == "Case.1"], ci = T, 
	ci.type = "lower", conf.level = 0.9)$interval$limits["LCL"]
ci.lower.limit.95.Pct <- enorm(Arsenic[Source == "Case.1"], ci = T, 
	ci.type = "lower", conf.level = 0.95)$interval$limits["LCL"]
ci.upper.limit.90.Pct <- enorm(Arsenic[Source == "Case.1"], ci = T, 
	ci.type = "upper", conf.level = 0.9)$interval$limits["UCL"]
ci.upper.limit.95.Pct <- enorm(Arsenic[Source == "Case.1"], ci = T, 
	ci.type = "upper", conf.level = 0.95)$interval$limits["UCL"]
Case.1.limits <- c(ci.lower.limit.90.Pct, ci.lower.limit.95.Pct, 
	ci.upper.limit.90.Pct, ci.upper.limit.95.Pct)

# Case 2
ci.lower.limit.90.Pct <- enorm(Arsenic[Source == "Case.2"], ci = T, 
	ci.type = "lower", conf.level = 0.9)$interval$limits["LCL"]
ci.lower.limit.95.Pct <- enorm(Arsenic[Source == "Case.2"], ci = T, 
	ci.type = "lower", conf.level = 0.95)$interval$limits["LCL"]
ci.upper.limit.90.Pct <- enorm(Arsenic[Source == "Case.2"], ci = T, 
	ci.type = "upper", conf.level = 0.9)$interval$limits["UCL"]
ci.upper.limit.95.Pct <- enorm(Arsenic[Source == "Case.2"], ci = T, 
	ci.type = "upper", conf.level = 0.95)$interval$limits["UCL"]
Case.2.limits <- c(ci.lower.limit.90.Pct, ci.lower.limit.95.Pct, 
	ci.upper.limit.90.Pct, ci.upper.limit.95.Pct)

mat <- round(rbind(Case.1.limits, Case.2.limits), 1)
dimnames(mat) <- list(c("Case.1", "Case.2"), 
	c("90%.LCL", "95%.LCL", "90%.UCL", "95%.UCL"))
mat
#       90%.LCL 95%.LCL 90%.UCL 95%.UCL
#Case.1    60.4    56.5    78.7    82.7
#Case.2    73.7    69.7    92.1    96.2


rm(Arsenic, Source, stats, kappa, UPL, UTL, 
	ci.lower.limit.90.Pct, ci.lower.limit.95.Pct,
	ci.upper.limit.90.Pct, ci.upper.limit.95.Pct,
	Case.1.limits, Case.2.limits, mat)


