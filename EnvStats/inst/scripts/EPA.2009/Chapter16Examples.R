#####
# File:     Script - EPA09, Chapter 16 Examples.R
#
# Purpose:  Reproduce Examples in Chapter 16 of the EPA Guidance Document
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
# Updated:  September 15, 2012
#####

#NOTE:  Unless otherwise noted, differences between what is shown in the
#       EPA Guidance Document and results shown here are due to the 
#       EPA Guidance Document rounding intermediate results prior to 
#       the final result.

#######################################################################################

library(EnvStats)

# Figure 16-1, p. 16-2
#---------------------

df.vec <- c(1, 3, 7, 25)
lty.vec <- c(1, 5, 4, 6)
col.vec <- 1:4
n <- length(df.vec)

windows()
pdfPlot(dist = "t", param.list = list(df = df.vec[1]), 
	curve.fill = F, xaxt = "n", 
	left.tail.cutoff = pt(-5, df = df.vec[1]), 
	right.tail.cutoff = pt(-5, df = df.vec[1]),
	xlab = "t-value", xlim = c(-5, 6),
	ylim = c(0, 0.4), main = "")
axis(1, at = seq(-5, 5, by = 2.5))

for(i in 2:n) {
	pdfPlot(dist = "t", param.list = list(df = df.vec[i]), 
		left.tail.cutoff = pt(-5, df = df.vec[i]), 
		right.tail.cutoff = pt(-5, df = df.vec[i]),
		add = T, pdf.col = col.vec[i], pdf.lty = lty.vec[i])
}

legend(3, 0.4, rev(paste(format(df.vec), "df")), lty = rev(lty.vec), 
	col = rev(col.vec), bty = "n", lwd = 2)

mtext(expression(paste("Figure 16-1. Student's ",
	italic(t), "-Distribtuion", sep = "")), cex = 1.25,
	line = 2.5)
mtext("for Varying Degrees of Freedom", cex = 1.25, 
	line = 1)

rm(df.vec, lty.vec, col.vec, n, i)


#######################################################################################

# Example 16-1, pp. 16-6 to 16-7
#-------------------------------

Sulfate   <- EPA.09.Ex.16.1.sulfate.df$Sulfate.ppm
Well.type <- EPA.09.Ex.16.1.sulfate.df$Well.type
fit    <- lm(Sulfate ~ Well.type, 
	data = EPA.09.Ex.16.1.sulfate.df)
Residuals <- fit$residuals

# The function summaryStats will compute summary statistics and
# also the p-value for testing the equality of means between groups.

summaryStats(Sulfate.ppm ~ Well.type, data = EPA.09.Ex.16.1.sulfate.df, 
	digits = 2, p.value = TRUE, alternative = "greater",
	stats.in.rows = TRUE)
#
#                        Background Downgradient Combined  
#N                         8          6           14       
#Mean                    536.25     608.33       567.14    
#SD                       26.69      18.35        43.4     
#Median                  540        605          565       
#Min                     490        590          490       
#Max                     570        630          630       
#p.value.between.greater                           0.000053
#95%.LCL.between                                  49.39    
#95%.UCL.between                                        Inf
#NA's                      0          2            2       
#N.Total                   8          8           16   


# The function stripChart will display the raw data by group, 
# show means and SDs for each group, and compute the p-value for 
# testing the equality of means between groups.

windows()
stripChart(Sulfate.ppm ~ Well.type, data = EPA.09.Ex.16.1.sulfate.df, 
	ylab = "Sulfate (ppm)", location.scale.digits = 2,
	p.value = TRUE, alternative = "greater")



# Probability plot on residuals
windows()
qqPlot(Residuals, add.line = T, 
	ylab = "Sulfate residuals (ppm)", 
	main = "Figure 16-2. Probability Plot of \nCombined Sulfate Residuals")


# The R function t.test performs Student's t-test.
# By default, var.equal is set to FALSE, so you have to set it
# to TRUE for this example.
t.test(Sulfate[Well.type == "Downgradient"], 
	Sulfate[Well.type == "Background"],
	alternative = "greater", var.equal = TRUE)
#
#Results of Hypothesis Test
#--------------------------
#
#Null Hypothesis:                 difference in means = 0
#
#Alternative Hypothesis:          True difference in means is greater than 0
#
#Test Name:                        Two Sample t-test
#
#Estimated Parameter(s):          mean of x = 608.3333
#                                 mean of y = 536.2500
#
#Data:                            Sulfate[Well.type == "Downgradient"] and Sulfate[Well.type == "Background"]
#
#Test Statistic:                  t = 5.660985
#
#P-value:                         5.271532e-05
#
#95% Confidence Interval:         LCL = 49.38883
#                                 UCL =      Inf


rm(Sulfate, Well.type, fit, Residuals)


#####################################################################################################

# Example 16-2, pp. 16-9 to 16-10
#--------------------------------

longToWide(EPA.09.Ex.16.2.benzene.df, "Benzene.ppb", "Month", "Well.type")
#    Background Downgradient
#Jan        0.5          0.5
#Feb        0.8          0.7
#Mar        1.6          4.6
#Apr        1.8          2.0
#May        1.1         16.7
#Jun       16.1         12.5
#Jul        1.6         26.3
#Aug        0.6        186.0


#NOTE: set var.equal=FALSE in the call to summaryStats
summaryStats(Benzene.ppb ~ Well.type, data = EPA.09.Ex.16.2.benzene.df,
	digits = 2, p.value = TRUE, test.arg.list = list(var.equal = FALSE),
	group.p.value.type = "between", alternative = "greater",
	stats.in.rows = TRUE)
#                              Background Downgradient Combined
#N                               8          8           16     
#Mean                            3.01      31.16        17.09  
#SD                              5.31      63.22        45.71  
#Median                          1.35       8.55         1.7   
#Min                             0.5        0.5          0.5   
#Max                            16.1      186          186     
#Welch.p.value.between.greater                           0.12  
#95%.LCL.between                                       -14.26  
#95%.UCL.between                                          Inf


Well.type <- EPA.09.Ex.16.2.benzene.df$Well.type
Benzene   <- EPA.09.Ex.16.2.benzene.df$Benzene.ppb

t.test(Benzene[Well.type == "Downgradient"], 
	Benzene[Well.type == "Background"],
	alternative = "greater", var.equal = FALSE)
#Results of Hypothesis Test
#--------------------------
#
#Null Hypothesis:                 difference in means = 0
#
#Alternative Hypothesis:          True difference in means is greater than 0
#
#Test Name:                       Welch Two Sample t-test
#
#Estimated Parameter(s):          mean of x = 31.1625
#                                 mean of y =  3.0125
#
#Data:                            Benzene[Well.type == "Downgradient"] and Benzene[Well.type == "Background"]
#
#Test Statistic:                  t = 1.254938
#
#Test Statistic Parameter:        df = 7.09878
#
#P-value:                         0.1246182
#
#95% Confidence Interval:         LCL = -14.25922
#                                 UCL =       Inf

set.seed(25)
windows()
stripChart(Benzene.ppb ~ Well.type, data = EPA.09.Ex.16.2.benzene.df,
	method = "jitter", ylim = c(-25, 200), 
	ylab = "Benzene (ppb)",
	p.value = TRUE, alternative = "greater", 
	test.arg.list = list(var.equal = FALSE))


rm(Well.type, Benzene)


#####################################################################################################

# Example 16-3, pp. 16-11 to 16-14
#---------------------------------

new.EPA.df <- EPA.09.Ex.16.2.benzene.df
new.EPA.df$Benzene.ppb <- log(new.EPA.df$Benzene.ppb)
names(new.EPA.df)[3] <- "log.Benzene.ppb"

table.df <- data.frame(
	longToWide(EPA.09.Ex.16.2.benzene.df, "Benzene.ppb", "Month", "Well.type"),
	longToWide(new.EPA.df, "log.Benzene.ppb", "Month", "Well.type")
	)
names(table.df)[3:4] <- c("log.Background", "log.Downgradient")
round(table.df, 3)
#    Background Downgradient log.Background log.Downgradient
#Jan        0.5          0.5         -0.693           -0.693
#Feb        0.8          0.7         -0.223           -0.357
#Mar        1.6          4.6          0.470            1.526
#Apr        1.8          2.0          0.588            0.693
#May        1.1         16.7          0.095            2.815
#Jun       16.1         12.5          2.779            2.526
#Jul        1.6         26.3          0.470            3.270
#Aug        0.6        186.0         -0.511            5.226


summaryStats(Benzene.ppb ~ Well.type, data = EPA.09.Ex.16.2.benzene.df)
#             N    Mean      SD Median Min   Max
#Background   8  3.0125  5.3108   1.35 0.5  16.1
#Downgradient 8 31.1625 63.2229   8.55 0.5 186.0


summaryStats(log(Benzene.ppb) ~ Well.type, data = EPA.09.Ex.16.2.benzene.df)
#             N   Mean     SD Median     Min    Max
#Background   8 0.3719 1.0825 0.2827 -0.6931 2.7788
#Downgradient 8 1.8757 1.9847 2.0259 -0.6931 5.2257
 

#Step 1 - Test raw data for normality
#------------------------------------

gofGroupTest(Benzene.ppb ~ Well.type, data = EPA.09.Ex.16.2.benzene.df)

#Results of Group Goodness-of-Fit Test
-------------------------------------
#
#Test Method:                     Wilk-Shapiro GOF (Normal Scores)
#
#Hypothesized Distribution:       Normal
#
#Data:                            Benzene.ppb
#
#Grouping Variable:               Well.type
#
#Data Source:                     EPA.09.Ex.16.2.benzene.df
#
#Number of Groups:                2
#
#Sample Sizes:                    Background   = 8
#                                 Downgradient = 8
#
#Test Statistic:                  z (G) = -5.801397
#
#P-values for
#Individual Tests:                Background   = 1.192816e-05
#                                 Downgradient = 3.459439e-05
#
#P-value for
#Group Test:                      3.288234e-09
#
#Alternative Hypothesis:          At least one group
#                                 does not come from a
#                                 Normal Distribution.



#Step 2 - Compute summary statistics and 
#         test log-transformed data for normality
#------------------------------------------------

summaryStats(log(Benzene.ppb) ~ Well.type,  
	data = EPA.09.Ex.16.2.benzene.df, 
	digits = 4, p.value = TRUE, 
	test.arg.list = list(var.equal = FALSE),
	alternative = "greater", stats.in.rows = TRUE)

#                              Background Downgradient Combined
#N                              8          8           16      
#Mean                           0.3719     1.8757       1.1238 
#SD                             1.0825     1.9847       1.7287 
#Median                         0.2827     2.0259       0.5289 
#Min                           -0.6931    -0.6931      -0.6931 
#Max                            2.7788     5.2257       5.2257 
#Welch.p.value.between.greater                          0.044  
#95%.LCL.between                                        0.0663 
#95%.UCL.between                                           Inf  


gofGroupTest(log(Benzene.ppb) ~ Well.type, data = EPA.09.Ex.16.2.benzene.df)

#Results of Group Goodness-of-Fit Test
#-------------------------------------
#
#Test Method:                     Wilk-Shapiro GOF (Normal Scores)
#
#Hypothesized Distribution:       Normal
#
#Data:                            log(Benzene.ppb)
#
#Grouping Variable:               Well.type
#
#Data Source:                     EPA.09.Ex.16.2.benzene.df
#
#Number of Groups:                2
#
#Sample Sizes:                    Background   = 8
#                                 Downgradient = 8
#
#Test Statistic:                  z (G) = -0.470135
#
#P-values for
#Individual Tests:                Background   = 0.04471158
#                                 Downgradient = 0.84933304
#
#P-value for
#Group Test:                      0.3191293
#
#Alternative Hypothesis:          At least one group
#                                 does not come from a
#                                 Normal Distribution.



# Steps 3-5 - Compute Welch t-test
#---------------------------------

Month       <- EPA.09.Ex.16.2.benzene.df$Month
Well.type   <- EPA.09.Ex.16.2.benzene.df$Well.type
Benzene     <- EPA.09.Ex.16.2.benzene.df$Benzene.ppb
log.Benzene <- log(Benzene)

t.test(log.Benzene[Well.type == "Downgradient"], 
	log.Benzene[Well.type == "Background"],
	alternative = "greater", var.equal = FALSE)
#
#Results of Hypothesis Test
#--------------------------
#
#Null Hypothesis:                 difference in means = 0
#
#Alternative Hypothesis:          True difference in means is greater than 0
#
#Test Name:                       Welch Two Sample t-test
#
#Estimated Parameter(s):          mean of x = 1.8757293
#                                 mean of y = 0.3718509
#
#Data:                            log.Benzene[Well.type == "Downgradient"] and log.Benzene[Well.type == "Background"]
#
#Test Statistic:                  t = 1.881485
#
#P-value:                         0.04352518
#
#95% Confidence Interval:         LCL = 0.06631003
#                                 UCL =        Inf


windows()
stripChart(log(Benzene.ppb) ~ Well.type, data = EPA.09.Ex.16.2.benzene.df,
	ylab = "log [ Benzene (ppb) ]",
	p.value = TRUE, alternative = "greater", 
	test.arg.list = list(var.equal = FALSE))


#Figure 16-3
#-----------
windows()
plot(as.numeric(Month), Benzene, type = "n", xaxt = "n",
	xlab = "Month", ylab = "Benzene (ppb)",
	main = "Figure 16-3. Benzene Time Series Plot")
axis(1, at = sort(unique(as.numeric(Month))), labels = levels(Month))

points(as.numeric(Month)[Well.type == "Background"], 
	Benzene[Well.type == "Background"],
	type = "b", pch = 3, lty = 2)
points(as.numeric(Month)[Well.type == "Downgradient"], 
	Benzene[Well.type == "Downgradient"],
	type = "b", pch = 15)

legend(1, 150, c("Downgradient", "Background"), 
	pch = c(15, 3), lty = 1:2)

rm(new.EPA.df, table.df, Month, Well.type, Benzene, log.Benzene)

#######################################################################################################

# Example 16-4, pp. 16-18 to 16-20
#---------------------------------

# NOTE: From the help file for wilcox.test:
#      "The literature is not unanimous about the definitions of the 
#       Wilcoxon rank sum and Mann-Whitney tests. The two most common 
#       definitions correspond to the sum of the ranks of the first sample 
#       with the minimum value subtracted or not: R subtracts and S-PLUS does not, 
#       giving a value which is larger by m(m+1)/2 for a first sample of size m. 
#       (It seems Wilcoxon's original paper used the unadjusted sum of the ranks 
#       but subsequent tables subtracted the minimum.) "


# Step 1
#-------

longToWide(EPA.09.Ex.16.4.copper.df, "Copper.ppb", "Month", "Well", paste.row.name = TRUE)
#        Well.1 Well.2 Well.3
#Month.1    4.2    5.2    9.4
#Month.2    5.8    6.4   10.1
#Month.3   11.3   11.3   14.5
#Month.4    7.0   11.5   16.1
#Month.5    7.0   10.1   21.5
#Month.6    8.2    9.7   17.6

new.EPA.df <- EPA.09.Ex.16.4.copper.df
new.EPA.df$Copper.ppb <- rank(new.EPA.df$Copper.ppb)
names(new.EPA.df)[4] <- "rank.Copper.ppb"
rank.df <- longToWide(new.EPA.df, "rank.Copper.ppb", "Month", "Well", 
	paste.row.name = TRUE)
rank.df
#        Well.1 Well.2 Well.3
#Month.1    1.0    2.0    8.0
#Month.2    3.0    4.0   10.5
#Month.3   12.5   12.5   15.0
#Month.4    5.5   14.0   16.0
#Month.5    5.5   10.5   18.0
#Month.6    7.0    9.0   17.0


# Step 2
#-------
colSums(rank.df)["Well.3"]
#Well.3 
#  84.5


# Step 3
#-------
Cu        <- EPA.09.Ex.16.4.copper.df$Copper.ppb
Well.type <- EPA.09.Ex.16.4.copper.df$Well.type

wilcox.test(Cu[Well.type == "Compliance"], Cu[Well.type == "Background"],
	alternative = "greater")

#Results of Hypothesis Test
#--------------------------
#
#Null Hypothesis:                 location shift = 0
#
#Alternative Hypothesis:          True location shift is greater than 0
#
#Test Name:                       Wilcoxon rank sum test with continuity correction
#
#Data:                            Cu[Well.type == "Compliance"] and Cu[Well.type == "Background"]
#
#Test Statistic:                  W = 63.5
#
#P-value:                         0.005659303
#
#Warning message:
#In wilcox.test.default(Cu[Well.type == "Compliance"], Cu[Well.type ==  :
#  cannot compute exact p-value with ties


rm(new.EPA.df, rank.df, Cu, Well.type)


#######################################################################################################

# Example 16-5, pp. 16-22 to 16-23
#---------------------------------

# NOTE: The EnvStats for R function 
#       twoSampleLinearRankTestCensored()
#       contains the Tarone-Ware test, along with several other
#       similar tests.

PCE <- EPA.09.Ex.16.5.PCE.df$PCE.ppb
Well.type <- EPA.09.Ex.16.5.PCE.df$Well.type
Censored <- EPA.09.Ex.16.5.PCE.df$Censored


# Tarone-Ware test with hypergeometric variance
#----------------------------------------------
twoSampleLinearRankTestCensored(
	x = PCE[Well.type == "Compliance"], 
	x.censored = Censored[Well.type == "Compliance"], 
	y = PCE[Well.type == "Background"], 
	y.censored = Censored[Well.type == "Background"], 
	test = "tarone.ware", 
	alternative = "greater", 
	variance = "hypergeometric") 

#Results of Hypothesis Test
#Based on Censored Data
#--------------------------
#
#Null Hypothesis:                 Fy(t) = Fx(t)
#
#Alternative Hypothesis:          Fy(t) > Fx(t) for at least one t
#
#Test Name:                       Two-Sample Linear Rank Test:
#                                 Tarone-Ware Test
#                                 with Hypergeometric Variance
#
#Censoring Side:                  left
#
#Data:                            x = PCE[Well.type == "Compliance"]
#                                 y = PCE[Well.type == "Background"]
#
#Censoring Variable:              x = Censored[Well.type == "Compliance"]
#                                 y = Censored[Well.type == "Background"]
#
#Sample Sizes:                    nx = 8
#                                 ny = 6
#
#Percent Censored:                x = 12.5%
#                                 y = 50.0%
#
#Test Statistics:                 nu     =  8.458912
#                                 var.nu = 20.912407
#                                 z      =  1.849748
#
#P-value:                         0.03217495



# Logrank test with hypergeometric variance
#------------------------------------------
twoSampleLinearRankTestCensored(
	x = PCE[Well.type == "Compliance"], 
	x.censored = Censored[Well.type == "Compliance"], 
	y = PCE[Well.type == "Background"], 
	y.censored = Censored[Well.type == "Background"], 
	test = "logrank", 
	alternative = "greater", 
	variance = "hypergeometric") 

#Results of Hypothesis Test
#Based on Censored Data
#--------------------------
#
#Null Hypothesis:                 Fy(t) = Fx(t)
#
#Alternative Hypothesis:          Fy(t) > Fx(t) for at least one t
#
#Test Name:                       Two-Sample Linear Rank Test:
#                                 Logrank Test
#                                 with Hypergeometric Variance
#
#Censoring Side:                  left
#
#Data:                            x = PCE[Well.type == "Compliance"]
#                                 y = PCE[Well.type == "Background"]
#
#Censoring Variable:              x = Censored[Well.type == "Compliance"]
#                                 y = Censored[Well.type == "Background"]
#
#Sample Sizes:                    nx = 8
#                                 ny = 6
#
#Percent Censored:                x = 12.5%
#                                 y = 50.0%
#
#Test Statistics:                 nu     = 2.830406
#                                 var.nu = 2.176723
#                                 z      = 1.918435
#
#P-value:                         0.02752793




# Peto-Peto test with hypergeometric variance
#--------------------------------------------
twoSampleLinearRankTestCensored(
	x = PCE[Well.type == "Compliance"], 
	x.censored = Censored[Well.type == "Compliance"], 
	y = PCE[Well.type == "Background"], 
	y.censored = Censored[Well.type == "Background"], 
	test = "peto.peto", 
	alternative = "greater", 
	variance = "hypergeometric")

#Results of Hypothesis Test
#Based on Censored Data
#--------------------------
#
#Null Hypothesis:                 Fy(t) = Fx(t)
#
#Alternative Hypothesis:          Fy(t) > Fx(t) for at least one t
#
#Test Name:                       Two-Sample Linear Rank Test:
#                                 Peto-Peto Test Using
#                                 Prentice Survival Estimator
#                                 with Hypergeometric Variance
#
#Censoring Side:                  left
#
#Data:                            x = PCE[Well.type == "Compliance"]
#                                 y = PCE[Well.type == "Background"]
#
#Censoring Variable:              x = Censored[Well.type == "Compliance"]
#                                 y = Censored[Well.type == "Background"]
#
#Sample Sizes:                    nx = 8
#                                 ny = 6
#
#Percent Censored:                x = 12.5%
#                                 y = 50.0%
#
#Test Statistics:                 nu     = 1.888889
#                                 var.nu = 1.028642
#                                 z      = 1.862406
#
#P-value:                         0.03127296


 
# Gehan test with hypergeometric variance
#----------------------------------------
twoSampleLinearRankTestCensored(
	x = PCE[Well.type == "Compliance"], 
	x.censored = Censored[Well.type == "Compliance"], 
	y = PCE[Well.type == "Background"], 
	y.censored = Censored[Well.type == "Background"], 
	test = "gehan", 
	alternative = "greater", 
	variance = "hypergeometric")

#Results of Hypothesis Test
#Based on Censored Data
#--------------------------
#
#Null Hypothesis:                 Fy(t) = Fx(t)
#
#Alternative Hypothesis:          Fy(t) > Fx(t) for at least one t
#
#Test Name:                       Two-Sample Linear Rank Test:
#                                 Gehan's Test
#                                 with Hypergeometric Variance
#
#Censoring Side:                  left
#
#Data:                            x = PCE[Well.type == "Compliance"]
#                                 y = PCE[Well.type == "Background"]
#
#Censoring Variable:              x = Censored[Well.type == "Compliance"]
#                                 y = Censored[Well.type == "Background"]
#
#Sample Sizes:                    nx = 8
#                                 ny = 6
#
#Percent Censored:                x = 12.5%
#                                 y = 50.0%
#
#Test Statistics:                 nu     =  27.000000
#                                 var.nu = 227.000000
#                                 z      =   1.792053
#
#P-value:                         0.03656224


rm(PCE, Well.type, Censored)

############################################################################