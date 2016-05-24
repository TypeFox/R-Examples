#####
# File:     Script - EPA09, Chapter 19 Examples.R
#
# Purpose:  Reproduce Examples in Chapter 19 of the EPA Guidance Document
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

# Example 19-1, pp. 19-17 to 19-18
#---------------------------------

names(EPA.09.Ex.19.1.sulfate.df)
#[1] "Well"                 "Month"                "Day"                 
#[4] "Year"                 "Date"                 "Sulfate.mg.per.l"    
#[7] "log.Sulfate.mg.per.l"

EPA.09.Ex.19.1.sulfate.df[, c("Well", "Date", "Sulfate.mg.per.l", "log.Sulfate.mg.per.l")]
#    Well       Date Sulfate.mg.per.l log.Sulfate.mg.per.l
#1  GW-01 1999-07-08             63.0             4.143135
#2  GW-01 1999-09-12             51.0             3.931826
#3  GW-01 1999-10-16             60.0             4.094345
#4  GW-01 1999-11-02             86.0             4.454347
#5  GW-04 1999-07-09            104.0             4.644391
#6  GW-04 1999-09-14            102.0             4.624973
#7  GW-04 1999-10-12             84.0             4.430817
#8  GW-04 1999-11-15             72.0             4.276666
#9  GW-08 1997-10-12             31.0             3.433987
#10 GW-08 1997-11-16             84.0             4.430817
#11 GW-08 1998-01-28             65.0             4.174387
#12 GW-08 1999-04-20             41.0             3.713572
#13 GW-08 2002-06-04             51.8             3.947390
#14 GW-08 2002-09-16             57.5             4.051785
#15 GW-08 2002-12-02             66.8             4.201703
#16 GW-08 2003-03-24             87.1             4.467057
#17 GW-09 1997-10-16             59.0             4.077537
#18 GW-09 1998-01-28             85.0             4.442651
#19 GW-09 1998-04-12             75.0             4.317488
#20 GW-09 1998-07-12             99.0             4.595120
#21 GW-09 2000-01-30             75.8             4.328098
#22 GW-09 2000-04-24             82.5             4.412798
#23 GW-09 2000-10-24             85.5             4.448516
#24 GW-09 2002-12-01            188.0             5.236442
#25 GW-09 2003-03-24            150.0             5.010635
 

# Summary statistics for raw data
#--------------------------------

summaryStats(Sulfate.mg.per.l ~ 1, 
	data = EPA.09.Ex.19.1.sulfate.df,
	digits = 2)
#                  N  Mean    SD Median Min Max
#Sulfate.mg.per.l 25 80.24 32.83   75.8  31 188


# Step 1:  Shapiro-Wilk Test on Background Data
#----------------------------------------------

gof.list <- gofTest(Sulfate.mg.per.l ~ 1, 
	data = EPA.09.Ex.19.1.sulfate.df)
gof.list
#
#Results of Goodness-of-Fit Test
#-------------------------------
#
#Test Method:                     Shapiro-Wilk GOF
#
#Hypothesized Distribution:       Normal
#
#Estimated Parameter(s):          mean = 80.24000
#                                 sd   = 32.82773
#
#Estimation Method:               mvue
#
#Data:                            Sulfate.mg.per.l
#
#Data Source:                     EPA.09.Ex.19.1.sulfate.df
#
#Sample Size:                     25
#
#Test Statistic:                  W = 0.8523182
#
#Test Statistic Parameter:        n = 25
#
#P-value:                         0.001949635
#
#Alternative Hypothesis:          True cdf does not equal the
#                                 Normal Distribution.

windows()
plot(gof.list)



# Try Lognormal transformation
gof.list.lnorm <- gofTest(Sulfate.mg.per.l ~ 1, 
	data = EPA.09.Ex.19.1.sulfate.df, 
	dist = "lnorm")
gof.list.lnorm
#
#Results of Goodness-of-Fit Test
#-------------------------------
#
#Test Method:                     Shapiro-Wilk GOF
#
#Hypothesized Distribution:       Lognormal
#
#Estimated Parameter(s):          meanlog = 4.3156194
#                                 sdlog   = 0.3756697
#
#Estimation Method:               mvue
#
#Data:                            Sulfate.mg.per.l
#
#Data Source:                     EPA.09.Ex.19.1.sulfate.df
#
#Sample Size:                     25
#
#Test Statistic:                  W = 0.9654866
#
#Test Statistic Parameter:        n = 25
#
#P-value:                         0.5340842
#
#Alternative Hypothesis:          True cdf does not equal the
#                                 Lognormal Distribution.

windows()
plot(gof.list.lnorm)




# Steps 2 and 3 - Determine Power for 1-of-2 and 1-of-3 retest designs assuming: 
#                 25 Background Measures, 
#                 10 Constituents, 
#                 50 Wells, 
#                  2 Evaluations per year, and
#                 Site Wide False Positive Rate (SWFPR) of 10%
#------------------------------------------------------------------------------

n  <- 25
nc <- 10
nw <- 50

# Set r = Number of Evaluations per year = 2
r <- 2


# Set Individual Test Type I Error to 
# 1 - (1 - SWFPR)^(1 / (Number of Constituents * Number of Wells))
#
# which translates to setting the confidence limit to 
# (1 - SWFPR)^(1 / (Number of Constituents * Number of Wells))

conf.level <- (1 - 0.1)^(1 / (nc * nw))
conf.level
#[1] 0.9997893


windows()
par(mar = c(5,4,6,1) + 0.1)
plotPredIntNormSimultaneousTestPowerCurve(n = n, k = 1, m = 3, 
	r = r, rule="k.of.m", pi.type = "upper", 
	conf.level = conf.level,  
	xlab = "SD Units Above Background",
	main = "")
plotPredIntNormSimultaneousTestPowerCurve(n = n, k = 1, m = 2, 
	r = r, rule="k.of.m", pi.type = "upper", 
	conf.level = conf.level,  
	add = T, plot.col = 2, plot.lty = 2)
legend(0, 1, c("1-of-3", "1-of-2"), col = 1:2, lty = 1:2, lwd = 2, bty="n")
mtext(paste('Power of 1-of-2 and 1-of-3 Retest Designs to Detect An Increase',
	"\nAt", nw, "Wells for SWFPR = 10%"), line = 3, cex = 1.25)
mtext(paste("(Adjusted for", nc, "Consituents and", r, "Evaluations per Year)"), line = 1)


# Multiplier for 1-of-2 Retest Design
predIntNormSimultaneousK(n = n, k = 1, m = 2,  
	r = r, rule = "k.of.m", pi.type = "upper", 
	conf.level = conf.level) 
#[1] 2.773504


# Multiplier for 1-of-3 Retest Design
predIntNormSimultaneousK(n = n, k = 1, m = 3,  
	r = r, rule = "k.of.m", pi.type = "upper", 
	conf.level = conf.level) 
#[1] 2.014364



# Step 4 - Compute upper prediction limit using 1-of-3 design
#          and assuming observations are lognormally distributed
#---------------------------------------------------------------

Sulfate <- EPA.09.Ex.19.1.sulfate.df$Sulfate.mg.per.l

pred.int.list <- 
	predIntLnormSimultaneous(x = Sulfate, 
		k = 1, m = 3, r = r, 
		rule = "k.of.m", pi.type = "upper", 
		conf.level = conf.level)

pred.int.list
#
#Results of Distribution Parameter Estimation
#--------------------------------------------
#
#Assumed Distribution:            Lognormal
#
#Estimated Parameter(s):          meanlog = 4.3156194
#                                 sdlog   = 0.3756697
#
#Estimation Method:               mvue
#
#Data:                            Sulfate
#
#Sample Size:                     25
#
#Prediction Interval Method:      exact 
#
#Prediction Interval Type:        upper
#
#Confidence Level:                99.97893%
#
#Minimum Number of
#Future Observations
#Interval Should Contain
#(per Sampling Occasion):         1
#
#Total Number of
#Future Observations
#(per Sampling Occasion):         3
#
#Number of Future
#Sampling Occasions:              100
#
#Prediction Interval:             LPL =   0.0000
#                                 UPL = 159.5496


names(pred.int.list)
#[1] "distribution" "sample.size"  "parameters"   "n.param.est" 
#[5] "method"       "data.name"    "bad.obs"      "interval"

pred.int.list$interval$limits["UPL"]
#     UPL 
#159.5496

# NOTE:  The upper prediction limit is on the scale of the original data,
#        *NOT* the log-transformed scale.

  
rm(Sulfate, gof.list, gof.list.lnorm, n, nc, nw, r, conf.level, pred.int.list)


#######################################################################################

# Example 19-2, pp. 19-18 to 19-20
#---------------------------------

# Raw data and summary statistics
#--------------------------------

head(EPA.09.Ex.19.2.chloride.df)
#   Well Chloride.mg.per.l
#1 GW-09              22.0
#2 GW-09              18.4
#3 GW-09              39.9
#4 GW-09              33.7
#5 GW-12              78.0
#6 GW-12              70.0

EPA.mat <- with(EPA.09.Ex.19.2.chloride.df,
	sapply(split(Chloride.mg.per.l, Well), I))
EPA.mat 
#     GW-09 GW-12 GW-13 GW-14 GW-15 GW-16 GW-24 GW-25 GW-26 GW-28
#[1,]  22.0  78.0  75.1  59.2  35.0  31.0  23.4  33.5  79.8  37.7
#[2,]  18.4  70.0  65.6  57.1  56.8  34.6  36.4  30.2  61.3  26.6
#[3,]  39.9  61.0  67.0  41.1  69.8  60.1  31.1  23.1  57.8  45.7
#[4,]  33.7  65.8  55.3  47.7  41.3  48.7  45.0  38.7  44.8  42.0

summaryStats(Chloride.mg.per.l ~ Well,
	data = EPA.09.Ex.19.2.chloride.df,  
	digits = 3, stats.in.rows = TRUE)
#       GW-09  GW-12  GW-13  GW-14  GW-15  GW-16  GW-24  GW-25  GW-26  GW-28 
#N       4      4      4      4      4      4      4      4      4      4    
#Mean   28.5   68.7   65.75  51.275 50.725 43.6   33.975 31.375 60.925 38    
#SD     10.021  7.208  8.128  8.427 15.672 13.392  9.083  6.533 14.447  8.273
#Median 27.85  67.9   66.3   52.4   49.05  41.65  33.75  31.85  59.55  39.85 
#Min    18.4   61     55.3   41.1   35     31     23.4   23.1   44.8   26.6  
#Max    39.9   78     75.1   59.2   69.8   60.1   45     38.7   79.8   45.7  
 
 

# Step 1 - Determine Power for Modified CA retest design assuming:
#           4 Background measures per Well, 
#           5 Constituents, 
#          10 Wells, 
#           1 Evaluation per year, and
#          Site Wide False Positive Rate (SWFPR) of 10%
#--------------------------------------------------------

n  <- 4
nc <- 5
nw <- 10

# Set r = Number of Evaluations per year = 1
r <- 1


# Set Individual Test Type I Error to 
# 1 - (1 - SWFPR)^(1 / (Number of Constituents * Number of Wells))
#
# which translates to setting the confidence limit to 
# (1 - SWFPR)^(1 / (Number of Constituents * Number of Wells))

conf.level <- (1 - 0.1)^(1 / (nc * nw))
conf.level
#[1] 0.997895


windows()
par(mar = c(5,4,6,1) + 0.1)
plotPredIntNormSimultaneousTestPowerCurve(n = n,  
	r = r, rule="Modified.CA", pi.type = "upper", 
	conf.level = conf.level,  
	ylim = c(0, 1),
	xlab = "SD Units Above Background",
	main = "")
mtext(paste('Power of Modified CA Retest Design to Detect An Increase',
	"\nAt", nw, "Wells for SWFPR = 10%"), line = 3, cex = 1.25)
mtext(paste("(Based on", n, "Background Observations per Well and\nAdjusted for", 
	nc, "Consituents and", r, "Evaluation per year)"), line = 0.5)


# Multiplier for Modified CA Retest Design
predIntNormSimultaneousK(n = n, r = r, rule = "Modified.CA", 
	pi.type = "upper", conf.level = conf.level) 
#[1] 4.356389



# Step 2 - Look at variability for each well
#-------------------------------------------
windows()
stripChart(Chloride.mg.per.l ~ Well, 
	data = EPA.09.Ex.19.2.chloride.df,
	ci.offset = 1, ylim = c(0, 100), 
	xlab = "Well", ylab = "Chloride (mg/l)",
	location.scale.text.cex = 0.6, 
	cex.axis = 0.85, n.text.cex = 0.85)
mtext("Chloride by Well with 95% CI's", line = 2.5, cex = 1.25, font = 2)



# Use Levene's test to test for difference in variances between wells
#--------------------------------------------------------------------
varGroupTest(Chloride.mg.per.l ~ Well, 
	data = EPA.09.Ex.19.2.chloride.df)
#
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
#Estimated Parameter(s):          GW-09 = 100.42000
#                                 GW-12 =  51.96000
#                                 GW-13 =  66.07000
#                                 GW-14 =  71.01583
#                                 GW-15 = 245.62250
#                                 GW-16 = 179.34000
#                                 GW-24 =  82.50917
#                                 GW-25 =  42.67583
#                                 GW-26 = 208.72917
#                                 GW-28 =  68.44667
#
#Data:                            Chloride.mg.per.l
#
#Grouping Variable:               Well
#
#Data Source:                     EPA.09.Ex.19.2.chloride.df
#
#Sample Sizes:                    GW-09 = 4
#                                 GW-12 = 4
#                                 GW-13 = 4
#                                 GW-14 = 4
#                                 GW-15 = 4
#                                 GW-16 = 4
#                                 GW-24 = 4
#                                 GW-25 = 4
#                                 GW-26 = 4
#                                 GW-28 = 4
#
#Test Statistic:                  F = 1.067300
#
#Test Statistic Parameters:       num df   =  9
#                                 denom df = 30
#
#P-value:                         0.4139585



# Step 3 - One-Way ANOVA to Compute Pooled SD Estimate
#-----------------------------------------------------
aov.fit <- aov(Chloride.mg.per.l ~ Well, 
	data = EPA.09.Ex.19.2.chloride.df)
anova.aov.fit <- anova(aov.fit)
anova.aov.fit
#            Df Sum Sq Mean Sq F value    Pr(>F)    
#Well         9 7585.3  842.81  7.5467 1.099e-05 ***
#Residuals   30 3350.4  111.68 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

Sp <- sqrt(anova.aov.fit["Residuals", "Mean Sq"])
df <- anova.aov.fit["Residuals", "Df"]
c(Sp, df)
#[1] 10.56782 30.00000


# Compare power of Modified CA retest design based on using 
# 4 background samples at each well
# resulting in an estimate of SD with 
# 3 degrees of freedom at each well 
#
# versus 
#
# pooling estimates of SD from all 10 wells in order to 
# use an estimate of SD based on 30 degrees of freedom

windows()
par(mar = c(5,4,6,1) + 0.1)
plotPredIntNormSimultaneousTestPowerCurve(n = n,  
	r = r, rule="Modified.CA", pi.type = "upper", 
	conf.level = conf.level,  
	ylim = c(0, 1),
	xlab = "SD Units Above Background",
	main = "")
plotPredIntNormSimultaneousTestPowerCurve(n = n,  
	df = df, r = r, rule="Modified.CA", 
	pi.type = "upper", conf.level = conf.level, 
	plot.col = 2, plot.lty = 2, add = T)
legend(0, 1, c("Pooled Well SD (30 df)", "Individual Well SD (3 df)"), 
	col = 2:1, lty = 2:1, lwd = 2, bty = 'n')
mtext(paste('Power of Modified CA Retest Design to Detect An Increase',
	"\nAt", nw, "Wells for SWFPR = 10%"), line = 3, cex = 1.25)
mtext(paste("(Based on", n, "Background Observations per Well and\nAdjusted for", 
	nc, "Consituents and", r, "Evaluation per year)"), line = 0.5)



# Step 4 - Compute K-multiplier based on using pooled estimate of SD
#-------------------------------------------------------------------

K <- predIntNormSimultaneousK(n = n, df = df, r = r, 
	rule = "Modified.CA", pi.type = "upper", 
	conf.level = conf.level) 
K
#[1] 1.980025


# Compute upper prediction limits for each well:
Means <- colMeans(EPA.mat)
Means
# GW-09  GW-12  GW-13  GW-14  GW-15  GW-16  GW-24  GW-25  GW-26  GW-28 
#28.500 68.700 65.750 51.275 50.725 43.600 33.975 31.375 60.925 38.000 
 
UPLs <- Means + K * Sp
round(UPLs, 1)
#GW-09 GW-12 GW-13 GW-14 GW-15 GW-16 GW-24 GW-25 GW-26 GW-28 
# 49.4  89.6  86.7  72.2  71.6  64.5  54.9  52.3  81.8  58.9



rm(EPA.mat, n, nc, nw, r, conf.level, aov.fit, anova.aov.fit, 
	Sp, df, K, Means, UPLs)


#####################################################################################################

# Example 19-3, pp. 19-23 to 19-24
#---------------------------------


# Determine Power for:
#   1-of-2, 1-of-3, 1-of-4, Modified CA designs based on single future observations
#   1-of-1, 1-of-2                      designs based on mean of 2 future observations
#   1-of-1                              design  based on mean of 3 future observations
#
# assuming: 
#                  25 Background Measures, 
#                  20 Constituents, 
#                 100 Wells, 
#                  2 Evaluations per year, and
#                 Site Wide False Positive Rate (SWFPR) of 10%
#------------------------------------------------------------------------------

n  <- 25
nc <- 20
nw <- 100

# Set r = Number of Evaluations per year = 2
r <- 2

# Set Individual Test Type I Error to 
# 1 - (1 - SWFPR)^(1 / (Number of Constituents * Number of Wells))
#
# which translates to setting the confidence limit to 
# (1 - SWFPR)^(1 / (Number of Constituents * Number of Wells))

conf.level <- (1 - 0.1)^(1 / (nc * nw))
conf.level
#[1] 0.9999473



# Power for designs based on SINGLE future observations
#------------------------------------------------------

windows()
par(mar = c(5,4,6,1) + 0.1)
plotPredIntNormSimultaneousTestPowerCurve(n = n, 
	k = 1, m = 4, r = r, rule="k.of.m", pi.type = "upper", 
	conf.level = conf.level,  
	xlab = "SD Units Above Background",
	main = "")

plotPredIntNormSimultaneousTestPowerCurve(n = n, 
	k = 1, m = 3, r = r, rule="k.of.m", pi.type = "upper", 
	conf.level = conf.level,  
	add = T, plot.col = 2, plot.lty = 2)

plotPredIntNormSimultaneousTestPowerCurve(n = n, 
	k = 1, m = 2, r = r, rule="k.of.m", pi.type = "upper", 
	conf.level = conf.level,  
	add = T, plot.col = 3, plot.lty = 3)

plotPredIntNormSimultaneousTestPowerCurve(n = n, 
	r = r, rule="Modified.CA", pi.type = "upper", 
	conf.level = conf.level,  
	add = T, plot.col = 4, plot.lty = 4)

legend(0, 1, c("1-of-4", "Modified CA", "1-of-3", "1-of-2"), 
	col = c(1, 4, 2, 3), lty = c(1, 4, 2, 3), lwd = 2, bty="n")
mtext(paste('Power of 1-of-2, 1-of-3, 1-of-4, and Modified CA Designs to',
	"\nDetect An Increase At", nw, "Wells for SWFPR = 10%"), 
	line = 3, cex = 1.25)
mtext(paste("(Based on", n, "Background Observations and\n", 
	"Adjusted for", nc, "Consituents and", r, "Evaluations per year)"), 
	line = 1)


# Power for designs based on the MEAN of future observations
#-----------------------------------------------------------

windows()
par(mar = c(5,4,6,1) + 0.1)
plotPredIntNormSimultaneousTestPowerCurve(n = n, 
	k = 1, m = 2, n.mean = 2, r = r, rule="k.of.m", 
	pi.type = "upper", conf.level = conf.level,  
	xlab = "SD Units Above Background",
	main = "")

plotPredIntNormSimultaneousTestPowerCurve(n = n, 
	k = 1, m = 1, n.mean = 2, r = r, rule="k.of.m", 
	pi.type = "upper", 
	conf.level = conf.level, 
	add = T, plot.col = 2, plot.lty = 2)

plotPredIntNormSimultaneousTestPowerCurve(n = n, 
	k = 1, m = 1, n.mean = 3, r = r, rule="k.of.m", 
	pi.type = "upper", conf.level = conf.level, 
	add = T, plot.col = 3, plot.lty = 3)

legend(0, 1, c("1-of-2, Order 2", "1-of-1, Order 3", "1-of-1, Order 2"), 
	col = c(1, 3, 2), lty = c(1, 3, 2), lwd = 2, 
	bty="n")
mtext(paste('Power of 1-of-1 and 1-of-2 Designs Based on Means to Detect',
	"\nAn Increase At", nw, "Wells for SWFPR = 10%"), 
	line = 3, cex = 1.25)
mtext(paste("(Based on", n, "Background Observations and\n", 
	"Adjusted for", nc, "Consituents and", r, "Evaluations per Year)"), 
	line = 1)



# Create data.frame containing K-multipliers and 
# Powers to detect an increase of 3 SDs above backgroud at all wells
#-------------------------------------------------------------------

rule.vec <- c(rep("k.of.m", 3), "Modified.CA", rep("k.of.m", 3))
m.vec <- c(2, 3, 4, 4, 1, 2, 1)
n.mean.vec <- c(rep(1, 4), 2, 2, 3)
n.scenarios <- length(rule.vec)
K.vec <- numeric(n.scenarios)
Power.vec <- numeric(n.scenarios)


for (i in 1:n.scenarios){
	K.vec[i] <- predIntNormSimultaneousK(
		n = n, k = 1, m = m.vec[i],  
		n.mean = n.mean.vec[i],
		r = r, rule = rule.vec[i], 
		pi.type = "upper", conf.level = conf.level)
	Power.vec[i] <- predIntNormSimultaneousTestPower(
		n = n, k = 1, m = m.vec[i], 
		n.mean = n.mean.vec[i], r = r, 
		rule = rule.vec[i], delta.over.sigma = 3, 
		pi.type = "upper", conf.level = conf.level)
}

data.frame(Rule = rule.vec, k = rep(1, n.scenarios), 
	m = m.vec, N.Mean = n.mean.vec, 
	K = round(K.vec, 2), 
	Power = round(Power.vec, 2), 
	Total.Samples = c(2, 3, 4, 4, 2, 4, 3))

#         Rule k m N.Mean    K Power Total.Samples
#1      k.of.m 1 2      1 3.16  0.39             2
#2      k.of.m 1 3      1 2.33  0.65             3
#3      k.of.m 1 4      1 1.83  0.81             4
#4 Modified.CA 1 4      1 2.57  0.71             4
#5      k.of.m 1 1      2 3.62  0.41             2
#6      k.of.m 1 2      2 2.33  0.85             4
#7      k.of.m 1 1      3 2.99  0.71             3


#NOTE:  The values of k and m have special meanings 
#       for the Modified CA rule.

rm(n, nc, nw, r, conf.level, rule.vec, m.vec, n.scenarios, 
	n.mean.vec, K.vec, Power.vec, i)


#######################################################################################

# Example 19-4, pp. 19-24 to 19-25
#---------------------------------

EPA.mat <- with(EPA.09.Ex.19.2.chloride.df,
	sapply(split(Chloride.mg.per.l, Well), I))
EPA.mat 
#     GW-09 GW-12 GW-13 GW-14 GW-15 GW-16 GW-24 GW-25 GW-26 GW-28
#[1,]  22.0  78.0  75.1  59.2  35.0  31.0  23.4  33.5  79.8  37.7
#[2,]  18.4  70.0  65.6  57.1  56.8  34.6  36.4  30.2  61.3  26.6
#[3,]  39.9  61.0  67.0  41.1  69.8  60.1  31.1  23.1  57.8  45.7
#[4,]  33.7  65.8  55.3  47.7  41.3  48.7  45.0  38.7  44.8  42.0


#  4 Background measures per well,
#  5 Constituents
# 10 Wells
#  1 Evaluation per year
n  <-  4
nc <-  5
nw <- 10
r  <-  1

# Set Individual Test Type I Error to 
# 1 - (1 - SWFPR)^(1 / (Number of Constituents * Number of Wells))
#
# which translates to setting the confidence limit to 
# (1 - SWFPR)^(1 / (Number of Constituents * Number of Wells))

conf.level <- (1 - 0.1)^(1 / (nc * nw))
conf.level
#[1] 0.997895


# Step 1 - One-Way ANOVA to Compute Pooled SD Estimate
#-----------------------------------------------------
aov.fit <- aov(Chloride.mg.per.l ~ Well, 
	data = EPA.09.Ex.19.2.chloride.df)
anova.aov.fit <- anova(aov.fit)
anova.aov.fit
#            Df Sum Sq Mean Sq F value    Pr(>F)    
#Well         9 7585.3  842.81  7.5467 1.099e-05 ***
#Residuals   30 3350.4  111.68 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

Sp <- sqrt(anova.aov.fit["Residuals", "Mean Sq"])
df <- anova.aov.fit["Residuals", "Df"]
c(Sp, df)
#[1] 10.56782 30.00000



# Step 2 - Compute K-multiplier based on 
#          using pooled estimate of SD and 
#          1-of-1, 1-of-2, and 1-of-3 plans using 
#          the mean of two future observations for
#          each evaluation
#-------------------------------------------------

m.vec <- 1:3
K.vec <- numeric(3)

for(i in 1:3) {
	K.vec[i] <- predIntNormSimultaneousK(
		n = n, df = df, k = 1, m = m.vec[i], 
		n.mean = 2, r = r, rule = "k.of.m", 
		pi.type = "upper", conf.level = conf.level)
} 

K.vec
#[1] 2.682613 1.875368 1.511350


# Compute upper prediction limits for each well:
Means <- colMeans(EPA.mat)
Means
# GW-09  GW-12  GW-13  GW-14  GW-15  GW-16  GW-24  GW-25  GW-26  GW-28 
#28.500 68.700 65.750 51.275 50.725 43.600 33.975 31.375 60.925 38.000
n.Wells <- length(Means) 
 
UPLs <- rep(Means, each = 3) + rep(K.vec, n.Wells) * Sp

data.frame(Well = names(UPLs), 
	Rule = rep(paste("1-of-", m.vec, sep=""), n.Wells),
	K = rep(round(K.vec, 2), n.Wells), 
	Mean = round(Means, 2), UPL = round(UPLs, 2))
#    Well   Rule    K  Mean   UPL
#1  GW-09 1-of-1 2.68 28.50 56.85
#2  GW-09 1-of-2 1.88 68.70 48.32
#3  GW-09 1-of-3 1.51 65.75 44.47
#4  GW-12 1-of-1 2.68 51.28 97.05
#5  GW-12 1-of-2 1.88 50.72 88.52
#6  GW-12 1-of-3 1.51 43.60 84.67
#7  GW-13 1-of-1 2.68 33.98 94.10
#8  GW-13 1-of-2 1.88 31.38 85.57
#9  GW-13 1-of-3 1.51 60.92 81.72
#10 GW-14 1-of-1 2.68 38.00 79.62
#11 GW-14 1-of-2 1.88 28.50 71.09
#12 GW-14 1-of-3 1.51 68.70 67.25
#13 GW-15 1-of-1 2.68 65.75 79.07
#14 GW-15 1-of-2 1.88 51.28 70.54
#15 GW-15 1-of-3 1.51 50.72 66.70
#16 GW-16 1-of-1 2.68 43.60 71.95
#17 GW-16 1-of-2 1.88 33.98 63.42
#18 GW-16 1-of-3 1.51 31.38 59.57
#19 GW-24 1-of-1 2.68 60.92 62.32
#20 GW-24 1-of-2 1.88 38.00 53.79
#21 GW-24 1-of-3 1.51 28.50 49.95
#22 GW-25 1-of-1 2.68 68.70 59.72
#23 GW-25 1-of-2 1.88 65.75 51.19
#24 GW-25 1-of-3 1.51 51.28 47.35
#25 GW-26 1-of-1 2.68 50.72 89.27
#26 GW-26 1-of-2 1.88 43.60 80.74
#27 GW-26 1-of-3 1.51 33.98 76.90
#28 GW-28 1-of-1 2.68 31.38 66.35
#29 GW-28 1-of-2 1.88 60.92 57.82
#30 GW-28 1-of-3 1.51 38.00 53.97
#Warning message:
#In data.frame(Well = names(UPLs), Rule = rep(paste("1-of-", m.vec,  :
#  row names were found from a short variable and have been discarded



# Step 3 - Figure 19-1
#---------------------

windows()
par(mar = c(5,4,6,1) + 0.1)
plotPredIntNormSimultaneousTestPowerCurve(
	n = n, df = df, k = 1, m = 3, n.mean = 2, r = r, 
	rule="k.of.m", pi.type = "upper", 
	conf.level = conf.level,  
	xlab = "SD Units Above Background",
	main = "")

plotPredIntNormSimultaneousTestPowerCurve(
	n = n, df = df, k = 1, m = 2, n.mean = 2, r = r, 
	rule="k.of.m", pi.type = "upper", 
	conf.level = conf.level, 
	add = T, plot.col = 2, plot.lty = 2)

plotPredIntNormSimultaneousTestPowerCurve(
	n = n, df =df, k = 1, m = 1, n.mean = 2, r = r, 
	rule="k.of.m", pi.type = "upper", 
	conf.level = conf.level, 
	add = T, plot.col = 3, plot.lty = 3)

legend(0, 1, c("1-of-3, Order 2", "1-of-2, Order 2", "1-of-1, Order 2"), 
	col = 1:3, lty = 1:3, lwd = 2, bty="n")
mtext(paste('Power of 1-of-1,2,3 Designs Based on Means to Detect',
	"\nAn Increase At", nw, "Wells for SWFPR = 10%"), 
	line = 3, cex = 1.25)
mtext(paste("(Based on", n, "Background Observations with", df, "DoF and\n", 
	"Adjusted for", nc, "Consituents and", r, "Evaluation per Year)"), 
	line = 1)



rm(EPA.mat, n, r, conf.level, aov.fit, anova.aov.fit, 
	Sp, df, K.vec, Means, n.Wells, UPLs)


############################################################################


# Example 19-5, pp. 19-33 to 19-35
#---------------------------------

#-------------------
# Show original data
#-------------------
head(EPA.09.Ex.19.5.mercury.df)
#  Event Well  Well.type Mercury.ppb.orig Mercury.ppb Censored
#1     1 BG-1 Background             0.21        0.21    FALSE
#2     2 BG-1 Background              <.2        0.20     TRUE
#3     3 BG-1 Background              <.2        0.20     TRUE
#4     4 BG-1 Background              <.2        0.20     TRUE
#5     5 BG-1 Background              <.2        0.20     TRUE
#6     6 BG-1 Background                           NA    FALSE


longToWide(EPA.09.Ex.19.5.mercury.df, "Mercury.ppb.orig", 
	"Event", "Well", paste.row.name = TRUE)
#        BG-1 BG-2 BG-3 BG-4 CW-1 CW-2
#Event.1 0.21  <.2  <.2  <.2 0.22 0.36
#Event.2  <.2  <.2 0.23 0.25  0.2 0.41
#Event.3  <.2  <.2  <.2 0.28  <.2 0.28
#Event.4  <.2 0.21 0.23  <.2 0.25 0.45
#Event.5  <.2  <.2 0.24  <.2 0.24 0.43
#Event.6                      <.2 0.54


#-------------------------------------
# Compute per-constituent Type I error
#-------------------------------------
# n  = 20 based on 4 Background Wells
# nw = 10 Compliance Wells
# nc =  5 Constituents
# ne =  1 Evaluation per year

n  <- 20
nw <- 10
nc <-  5
ne <-  1

# Set number of future sampling occasions r to 
# Number Compliance Wells x Number Evaluations per Year
r  <-  nw * ne

# Set Per-Constituent Type I Error to 
# 1 - (1 - SWFPR)^(1 / Number of Constituents)
#
# which translates to setting the confidence limit to 
# (1 - SWFPR)^(1 / Number of Constituents)

conf.level <- (1 - 0.1)^(1 / nc)
conf.level
#[1] 0.9791484

alpha <- 1 - conf.level
alpha
#[1] 0.02085164


#--------------------------------------------------
# Compute alpha levels for various resampling plans
#--------------------------------------------------

# Check alpha levels for the following retesting plans:
# 1) 1-of-2
# 2) 1-of-3
# 3) 1-of-4
# 4) Modified CA
# 5) 1-of-1 for median of 3 future values.  
#    This plan is equivalent to the 2-of-3 plan for single observations.
# 6) 1-of-2 for median of 3 future values.  

rule.vec <- c(rep("k.of.m", 3), "Modified.CA", rep("k.of.m", 2))
k.vec    <- rep(1, 6)
m.vec    <- c(2:4, 4, 1, 2)
n.median.vec <- c(rep(1, 4), rep(3, 2))

n.plans <- length(rule.vec)

# Alpha-levels based on using Maximum background observation
alpha.vec.Max <- numeric(n.plans)
for(i in 1:n.plans)
	alpha.vec.Max[i] <- 1 - 
		predIntNparSimultaneousConfLevel(
		n = n, n.median = n.median.vec[i], 
		k = k.vec[i], m = m.vec[i], r = r, 
		rule = rule.vec[i], pi.type="upper")

# Alpha-levels based on using 2nd largest background observation
alpha.vec.2nd <- numeric(n.plans)
for(i in 1:n.plans)
	alpha.vec.2nd[i] <- 1 - 
		predIntNparSimultaneousConfLevel(
		n = n, n.median = n.median.vec[i], 
		k = k.vec[i], m = m.vec[i], r = r, 
		rule = rule.vec[i], pi.type="upper",
		n.plus.one.minus.upl.rank = 2)

# Alpha-levels based on using 3rd largest background observation
alpha.vec.3rd <- numeric(n.plans)
for(i in 1:n.plans)
	alpha.vec.3rd[i] <- 1 - 
		predIntNparSimultaneousConfLevel(
		n = n, n.median = n.median.vec[i], 
		k = k.vec[i], m = m.vec[i], r = r, 
		rule = rule.vec[i], pi.type="upper",
		n.plus.one.minus.upl.rank = 3)

Well.type <- EPA.09.Ex.19.5.mercury.df$Well.type
Hg        <- EPA.09.Ex.19.5.mercury.df$Mercury.ppb
BG.sorted <- sort(Hg[Well.type == "Background"], decreasing = TRUE)
BG.sorted
# [1] 0.28 0.25 0.24 0.23 0.23 0.21 0.21 0.20 0.20 0.20 0.20 0.20 0.20 0.20
#[15] 0.20 0.20 0.20 0.20 0.20 0.20

Candidate.Plans.df <- data.frame(Rule = rep(rule.vec, 3), 
	Median.n = rep(n.median.vec, 3), 
	k = rep(k.vec, 3), m = rep(m.vec, 3),
	Order.Statistic = rep(c("Max", "2nd", "3rd"), each = n.plans), 
	Achieved.alpha = round(
		c(alpha.vec.Max, alpha.vec.2nd, alpha.vec.3rd), 4),
	BG.Limit = rep(BG.sorted[1:3], each = n.plans))

Candidate.Plans.df
#          Rule Median.n k m Order.Statistic Achieved.alpha BG.Limit
#1       k.of.m        1 1 2             Max         0.0395     0.28
#2       k.of.m        1 1 3             Max         0.0055     0.28
#3       k.of.m        1 1 4             Max         0.0009     0.28
#4  Modified.CA        1 1 4             Max         0.0140     0.28
#5       k.of.m        3 1 1             Max         0.0961     0.28
#6       k.of.m        3 1 2             Max         0.0060     0.28
#7       k.of.m        1 1 2             2nd         0.1118     0.25
#8       k.of.m        1 1 3             2nd         0.0213     0.25
#9       k.of.m        1 1 4             2nd         0.0046     0.25
#10 Modified.CA        1 1 4             2nd         0.0516     0.25
#11      k.of.m        3 1 1             2nd         0.2474     0.25
#12      k.of.m        3 1 2             2nd         0.0268     0.25
#13      k.of.m        1 1 2             3rd         0.2082     0.24
#14      k.of.m        1 1 3             3rd         0.0516     0.24
#15      k.of.m        1 1 4             3rd         0.0135     0.24
#16 Modified.CA        1 1 4             3rd         0.1170     0.24
#17      k.of.m        3 1 1             3rd         0.4166     0.24
#18      k.of.m        3 1 2             3rd         0.0709     0.24



#--------------------------------------------------------------------
# Eliminate plans that do not achieve the per-constituent alpha level
#--------------------------------------------------------------------
index <- Candidate.Plans.df$Achieved.alpha <= alpha
Candidate.Plans.df <- Candidate.Plans.df[index, ]

Candidate.Plans.df
#          Rule Median.n k m Order.Statistic Achieved.alpha BG.Limit
#2       k.of.m        1 1 3             Max         0.0055     0.28
#3       k.of.m        1 1 4             Max         0.0009     0.28
#4  Modified.CA        1 1 4             Max         0.0140     0.28
#6       k.of.m        3 1 2             Max         0.0060     0.28
#9       k.of.m        1 1 4             2nd         0.0046     0.25
#15      k.of.m        1 1 4             3rd         0.0135     0.24





#-------------------------------------------------
# Determine whether Well CW-1 is out of compliance
#-------------------------------------------------

EPA.09.Ex.19.5.mercury.df[Well == "CW-1", c("Well", "Event", "Mercury.ppb.orig")]
#   Well Event Mercury.ppb.orig
#25 CW-1     1             0.22
#26 CW-1     2              0.2
#27 CW-1     3              <.2
#28 CW-1     4             0.25
#29 CW-1     5             0.24
#30 CW-1     6              <.2


#        Rule Median.n k m Order.Statistic Achieved.alpha BG.Limit Result
# ----------- -------- - - --------------- -------------- -------- ------
#      k.of.m        1 1 3             Max         0.0055     0.28 Pass. First value of 0.22 < BG Limit.
#      k.of.m        1 1 4             Max         0.0009     0.28 Pass. First value of 0.22 < BG Limit.
# Modified.CA        1 1 4             Max         0.0140     0.28 Pass. First value of 0.22 < BG Limit.
#      k.of.m        3 1 2             Max         0.0060     0.28 Pass. First value of 0.22 < BG Limit.
#      k.of.m        1 1 4             2nd         0.0046     0.25 Pass. First value of 0.22 < BG Limit.
#      k.of.m        1 1 4             3rd         0.0135     0.24 Pass. First value of 0.22 < BG Limit.




#-------------------------------------------------
# Determine whether Well CW-2 is out of compliance
#-------------------------------------------------

EPA.09.Ex.19.5.mercury.df[Well == "CW-2", c("Well", "Event", "Mercury.ppb.orig")]
#   Well Event Mercury.ppb.orig
#31 CW-2     1             0.36
#32 CW-2     2             0.41
#33 CW-2     3             0.28
#34 CW-2     4             0.45
#35 CW-2     5             0.43
#36 CW-2     6             0.54


# Resampling is required for all plans because initial value of 0.36
# is greater than BG Limit, and median of first 3 values (also 0.36)
# is greater than BG Limit.

#        Rule Median.n k m Order.Statistic Achieved.alpha BG.Limit Result
# ----------- -------- - - --------------- -------------- -------- ------
#      k.of.m        1 1 3             Max         0.0055     0.28 Pass. Third value of 0.28 = BG Limit.
#      k.of.m        1 1 4             Max         0.0009     0.28 Pass. Third value of 0.28 = BG Limit.
# Modified.CA        1 1 4             Max         0.0140     0.28 Fail. Only 1 of observations 2-4 <= 0.28.
#      k.of.m        3 1 2             Max         0.0060     0.28 Fail. Both medians (0.36 and 0.45) > BG Limit.
#      k.of.m        1 1 4             2nd         0.0046     0.25 Fail. Observations 1-4 all > 0.25.
#      k.of.m        1 1 4             3rd         0.0135     0.24 Fail. Observations 1-4 all > 0.25.




#------------------------------------------------------------------------------
# Compute approximate power of each plan to detect contamination at just 1 well 
# assuming true underying distribution of Hg is Normal at all wells.
#------------------------------------------------------------------------------

windows()
plotPredIntNparSimultaneousTestPowerCurve(n = 20, 
	k = 1, m = 4, r = 10, rule = "k.of.m", 
	n.plus.one.minus.upl.rank = 3, 
	pi.type = "upper", r.shifted = 1, 
	method = "approx", range.delta.over.sigma = c(0, 5), main = "")

plotPredIntNparSimultaneousTestPowerCurve(n = 20, 
	n.median = 3, k = 1, m = 2, r = 10, rule = "k.of.m", 
	n.plus.one.minus.upl.rank = 1, 
	pi.type = "upper", r.shifted = 1, 
	method = "approx", range.delta.over.sigma = c(0, 5), 
	add = TRUE, plot.col = 2, plot.lty = 2)

plotPredIntNparSimultaneousTestPowerCurve(n = 20, 
	r = 10, rule = "Modified.CA", 
	n.plus.one.minus.upl.rank = 1, 
	pi.type = "upper", r.shifted = 1, 
	method = "approx", range.delta.over.sigma = c(0, 5), 
	add = TRUE, plot.col = 3, plot.lty = 3)

plotPredIntNparSimultaneousTestPowerCurve(n = 20, 
	k = 1, m = 4, r = 10, rule = "k.of.m", 
	n.plus.one.minus.upl.rank = 2, 
	pi.type = "upper", r.shifted = 1, 
	method = "approx", range.delta.over.sigma = c(0, 5), 
	add = TRUE, plot.col = 4, plot.lty = 4)

plotPredIntNparSimultaneousTestPowerCurve(n = 20, 
	k = 1, m = 3, r = 10, rule = "k.of.m", 
	n.plus.one.minus.upl.rank = 1, 
	pi.type = "upper", r.shifted = 1, 
	method = "approx", range.delta.over.sigma = c(0, 5), 
	add = TRUE, plot.col = 5, plot.lty = 5)

plotPredIntNparSimultaneousTestPowerCurve(n = 20, 
	k = 1, m = 4, r = 10, rule = "k.of.m", 
	n.plus.one.minus.upl.rank = 1, 
	pi.type = "upper", r.shifted = 1, 
	method = "approx", range.delta.over.sigma = c(0, 5), 
	add = TRUE, plot.col = 6, plot.lty = 6)

legend(0, 1, legend = c("1-of-4, 3rd", "1-of-2, Max, Median", "Mod CA", 
	"1-of-4, 2nd", "1-of-3, Max", "1-of-4, Max"), 
	lwd = 3, col = 1:6, lty = 1:6, bty = "n")

title(main = "Figure 19-2. Comparison of Full Power Curves")


rm(Well.type, Hg, n, nw, nc, ne, r, conf.level, alpha, 
	n.median.vec, rule.vec, k.vec, m.vec, n.plans, 
	alpha.vec.Max, alpha.vec.2nd, alpha.vec.3rd,
	BG.sorted, i, Candidate.Plans.df, index)
