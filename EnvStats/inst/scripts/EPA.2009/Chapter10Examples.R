#####
# File:     Script - EPA09, Chapter 10 Examples.R
#
# Purpose:  Reproduce Examples in Chapter 10 of the EPA Guidance Document
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
# Updated:  August 24, 2012
#####

#NOTE:  Unless otherwise noted, differences between what is shown in the
#       EPA Guidance Document and results shown here are due to the 
#       EPA Guidance Document rounding intermediate results prior to 
#       the final result.

#######################################################################################

library(EnvStats)

# Example 10-1, pp. 10-11 to 10-12
#---------------------------------

#NOTE: By default, EnvironmentalStats computes the skew using the 
#      Fisher method, whereas the example here uses the 
#      method of moments; thus the need to use the skew.list
#      argument in the call to summaryFull.

# Just Mean, SD, CV, and Skew:

summaryFull(Nickel.ppb ~ 1, data = EPA.09.Ex.10.1.nickel.df,
	stats = c("n", "mean", "sd", "cv", "skew"),
	skew.list = list(method="moment"), data.name = "Nickel")
#                         Nickel
#N                         20       
#Mean                     169.5     
#Skew                       1.843   
#Standard Deviation       259.7     
#Coefficient of Variation   1.532


# All summary stats:

summaryFull(Nickel.ppb ~ 1, data = EPA.09.Ex.10.1.nickel.df,
	skew.list = list(method="moment"), data.name = "Nickel")
#                             Nickel
#N                             20       
#Mean                         169.5     
#Median                        57.4     
#10% Trimmed Mean             113       
#Geometric Mean                50.33    
#Skew                           1.843   
#Kurtosis                       3.417   
#Min                            1       
#Max                          942       
#Range                        941       
#1st Quartile                  17.75    
#3rd Quartile                 178.8     
#Standard Deviation           259.7     
#Geometric Standard Deviation   6.058   
#Interquartile Range          161.1     
#Median Absolute Deviation     67.31    
#Coefficient of Variation       1.532


summaryFull(log(Nickel.ppb) ~ 1, data = EPA.09.Ex.10.1.nickel.df,
	stats = c("n", "mean", "sd", "cv", "skew"),
	skew.list = list(method="moment"), data.name = "log(Nickel)")
#                         log(Nickel)
#N                        20             
#Mean                      3.919         
#Skew                     -0.245         
#Standard Deviation        1.801         
#Coefficient of Variation  0.4597


################################################################


# Example 10-2, pp. 10-14 to 10-15
#---------------------------------

gofTest(EPA.09.Ex.10.1.nickel.df$Nickel, data.name.x = "Nickel")

#Results of Goodness-of-Fit Test
#-------------------------------
#
#Test Method:                     Shapiro-Wilk GOF
#
#Hypothesized Distribution:       Normal
#
#Estimated Parameter(s):          mean = 169.5250
#                                 sd   = 259.7175
#
#Estimation Method:               mvue
#
#Data:                            Nickel
#
#Sample Size:                     20
#
#Test Statistic:                  W = 0.6788888
#
#Test Statistic Parameter:        n = 20
#
#P-value:                         2.17927e-05
#
#Alternative Hypothesis:          True cdf does not equal the
#                                 Normal Distribution.


###################################################################


# Example 10-3, pp. 10-17 to 10-18
#---------------------------------

gofTest(EPA.09.Ex.10.1.nickel.df$Nickel, data.name.x = "Nickel", 
	test = "ppcc")

#Results of Goodness-of-Fit Test
#-------------------------------
#
#Test Method:                     PPCC GOF
#
#Hypothesized Distribution:       Normal
#
#Estimated Parameter(s):          mean = 169.5250
#                                 sd   = 259.7175
#
#Estimation Method:               mvue
#
#Data:                            Nickel
#
#Sample Size:                     20
#
#Test Statistic:                  r = 0.8199825
#
#Test Statistic Parameter:        n = 20
#
#P-value:                         5.753418e-05
#
#Alternative Hypothesis:          True cdf does not equal the
#                                 Normal Distribution.


######################################################################


# Example 10-4, pp. 10-20 to 10-21
#---------------------------------

#NOTE: The table shown on page 10-12 shows nickel
#      concentrations by Month and Year.  However, 
#      the EPA Errata Sheet states "Year" should be "Well".
#
#      The table on page 10-21 shows values that 
#      correspond to the log-transformed values shown in 
#      the table on page 10-12.  However, the table on page 10-21
#      is labeled with Well number for the columns (not Year number),
#      and Month is labeled as 1, 2, 3, 4, 5 instead of 
#      Jan, Mar, Jun, Aug, Oct.
#
#      In this script, we will use the data in the table 
#      shown on page 10-12, and do a group test by Well.
#
#      NOTE: The Shapiro-Wilk statistics by Well shown on 
#            page 10-20 are incorrect, as is the overall group statistic.
#            However, the EPA Errata Sheet shows correct statistics.
 

# Show Shapiro-Wilks statistics by Well

SW.stats <- sapply(split(EPA.09.Ex.10.1.nickel.df$Nickel.ppb, 
	EPA.09.Ex.10.1.nickel.df$Well), 
	function(x) gofTest(x)$statistic)
round(SW.stats, 4)
#Well.1.W Well.2.W Well.3.W Well.4.W 
#  0.7578   0.7397   0.7066   0.8150 

rm(SW.stats)


gofGroupTest(Nickel.ppb ~ Well, data = EPA.09.Ex.10.1.nickel.df)

#Results of Group Goodness-of-Fit Test
-------------------------------------
#
#Test Method:                     Wilk-Shapiro GOF (Normal Scores)
#
#Hypothesized Distribution:       Normal
#
#Data:                            Nickel.ppb
#
#Grouping Variable:               Well
#
#Data Source:                     EPA.09.Ex.10.1.nickel.df
#
#Number of Groups:                4
#
#Sample Sizes:                    Well.1 = 5
#                                 Well.2 = 5
#                                 Well.3 = 5
#                                 Well.4 = 5
#
#Test Statistic:                  z (G) = -3.658696
#
#P-values for
#Individual Tests:                Well.1 = 0.03510747
#                                 Well.2 = 0.02385344
#                                 Well.3 = 0.01120775
#                                 Well.4 = 0.10681461
#
#P-value for
#Group Test:                      0.0001267509
#
#Alternative Hypothesis:          At least one group
#                                 does not come from a
#                                 Normal Distribution.


###################################################################


# Example 10-5, pp. 10-21 to 10-24
#---------------------------------


summaryFull(log(Nickel.ppb) ~ 1, data = EPA.09.Ex.10.1.nickel.df, 
	stats = c("n", "mean", "sd", "cv", "skew"),
	skew.list = list(method="moment"))
#                         log(Nickel.ppb)
#N                        20             
#Mean                      3.919         
#Skew                     -0.245         
#Standard Deviation        1.801         
#Coefficient of Variation  0.4597


gofTest(EPA.09.Ex.10.1.nickel.df$Nickel.ppb, dist = "lnorm", 
	data.name.x = "Nickel")

#Results of Goodness-of-Fit Test
#-------------------------------
#
#Test Method:                     Shapiro-Wilk GOF
#
#Hypothesized Distribution:       Lognormal
#
#Estimated Parameter(s):          meanlog = 3.918529
#                                 sdlog   = 1.801404
#
#Estimation Method:               mvue
#
#Data:                            Nickel
#
#Sample Size:                     20
#
#Test Statistic:                  W = 0.978946
#
#Test Statistic Parameter:        n = 20
#
#P-value:                         0.9197735
#
#Alternative Hypothesis:          True cdf does not equal the
#                                 Lognormal Distribution.


gofTest(EPA.09.Ex.10.1.nickel.df$Nickel, test = "ppcc", 
	dist = "lnorm", data.name.x = "Nickel")

#Results of Goodness-of-Fit Test
#-------------------------------
#
#Test Method:                     PPCC GOF
#
#Hypothesized Distribution:       Lognormal
#
#Estimated Parameter(s):          meanlog = 3.918529
#                                 sdlog   = 1.801404
#
#Estimation Method:               mvue
#
#Data:                            Nickel
#
#Sample Size:                     20
#
#Test Statistic:                  r = 0.9912528
#
#Test Statistic Parameter:        n = 20
#
#P-value:                         0.9187852
#
#Alternative Hypothesis:          True cdf does not equal the
#                                 Lognormal Distribution.






