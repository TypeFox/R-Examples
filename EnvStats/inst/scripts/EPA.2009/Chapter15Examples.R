#####
# File:     Script - EPA09, Chapter 15 Examples.R
#
# Purpose:  Reproduce Examples in Chapter 15 of the EPA Guidance Document
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
# Updated:  2013/08/27
#####

#NOTE:  Unless otherwise noted, differences between what is shown in the
#       EPA Guidance Document and results shown here are due to the 
#       EPA Guidance Document rounding intermediate results prior to 
#       the final result.

#######################################################################################

library(EnvStats)

# Example 15-1, pp. 15-10 to 15-13
#---------------------------------

#------------------------
# Results for raw Mn data
#------------------------


# Figure 15-3
#------------

# NOTE:  Figure 15-3 displays quantiles based on plotting positions 
#        associated with the observed (detected) Manganese concentrations,
#        including the largest observation.  In order to display the 
#        plotting position associated with the largest observation 
#        when using the Kaplan-Meier method with left-censored data, USEPA (2009) 
#        modifies the K-M method for the largest value, because for the K-M method, 
#        the empirical cdf is equal to 1 for the largest value.  So in order to 
#        reproduce Figure 15-3, you have to set prob.method="kaplan-meier with max" 
#        in the call to qqplotCensored().  
#        Setting prob.method="kaplan-meier" will produce the exact same plot except 
#        there won't be any plotting position associated with the largest value.

#        Also, Figure 15-3 displays quantiles based on plotting positions 
#        associated with the censoring levels <2 and <5.
#        By default, qqPlotCensored does not display quantiles associated 
#        with the censoroing levels and this is the recommended method 
#        (see for example Helsel, 2011, p. 52).
#        Here, we will create two plots:
#          one without displaying quantiles associated with the censoring levels, 
#          one with    displaying quantiles associated with the censoring levels.
 

# Omit estimated quantiles associated with the censoring levels
#--------------------------------------------------------------
windows()
qqplot.list <- with(EPA.09.Ex.15.1.manganese.df, 
	qqPlotCensored(x = Manganese.ppb, censored = Censored, 
		prob.method = "modified kaplan-meier",
		ylim = c(0, max(Manganese.ppb)),
		add.line = TRUE, ylab = "manganese (ppb)", main = ""))
mtext(paste("Figure 15-3. Censored Probability Plot of",
	"Manganese Concentrations (Censoring Levels Not Shown)",
	sep = "\n"), cex = 1.15, font = 2, line = 1.75)

# Correlation between *DETECTED* observations and plotting positions:
#-------------------------------------------------------------------
cor(qqplot.list$x[!qqplot.list$Censored], qqplot.list$y[!qqplot.list$Censored])
#[1] 0.9132335


# Include estimated quantiles associated with the censoring levels
#-----------------------------------------------------------------
windows()
set.seed(566)
qqplot.list <- with(EPA.09.Ex.15.1.manganese.df, 
	qqPlotCensored(x = Manganese.ppb, censored = Censored, 
		prob.method = "modified kaplan-meier",
		ylim = c(0, max(Manganese.ppb)),
		include.cen = TRUE,  
		duplicate.points.method = "jitter",
		add.line = TRUE, ylab = "manganese (ppb)", main = ""))
mtext(paste("Figure 15-3. Censored Probability Plot of",
	"Manganese Concentrations", sep = "\n"), cex = 1.25, 
	font = 2, line = 1.75)

# Correlation between *UNIQUE* detected observations and censoring levels
#                     vs. plotting positions:
#------------------------------------------------------------------------
x <- qqplot.list$x
y <- qqplot.list$y

unique.y <- unique(y)
new.x <- x[match(unique.y, y)]

cor(new.x, unique.y)
#[1] 0.9021066



#------------------------------------
# Results for log-transformed Mn data
#------------------------------------


# Figure 15-4
#------------

# NOTE:  Figure 15-4 displays quantiles based on plotting positions 
#        associated with the observed (detected) Manganese concentrations,
#        including the largest observation.  In order to display the 
#        plotting position associated with the largest observation 
#        when using the Kaplan-Meier method with left-censored data, USEPA (2009) 
#        modifies the K-M method for the largest value, because for the K-M method, 
#        the empirical cdf is equal to 1 for the largest value.  So in order to 
#        reproduce Figure 15-4, you have to set prob.method="kaplan-meier with max" 
#        in the call to qqplotCensored().
#        Setting prob.method="kaplan-meier" will produce the exact same plot except 
#        there won't be any plotting position associated with the largest value.

#        Also, Figure 15-4 displays quantiles based on plotting positions 
#        associated with the censoring levels <2 and <5.
#        By default, qqPlotCensored does not display quantiles associated 
#        with the censoroing levels and this is the recommended method 
#        (see for example Helsel, 2011, p. 52).
#        Here, we will create two plots:
#          one without displaying quantiles associated with the censoring levels, 
#          one with    displaying quantiles associated with the censoring levels.

# Omit estimated quantiles associated with the censoring levels
#--------------------------------------------------------------
windows()
qqplot.list <- with(EPA.09.Ex.15.1.manganese.df, 
	qqPlotCensored(x = Manganese.ppb, censored = Censored, 
		prob.method = "modified kaplan-meier",
		dist = "lnorm", add.line = T, 
		ylab = "log(manganese) log(ppb)", main = ""))
mtext(paste("Figure 15-4. Censored Probability Plot of", 
	"Logged Manganese Sample (Censoring Levels Not Shown)", sep = "\n"), 
	cex = 1.15, font = 2, line = 1.75)

# Correlation between *DETECTED* observations and plotting positions:
#-------------------------------------------------------------------
cor(qqplot.list$x[!qqplot.list$Censored], qqplot.list$y[!qqplot.list$Censored])
#[1] 0.9931783


# Include estimated quantiles associated with the censoring levels
#-----------------------------------------------------------------
windows()
qqplot.list <- with(EPA.09.Ex.15.1.manganese.df,
	qqPlotCensored(x = Manganese.ppb, censored = Censored, 
		prob.method = "modified kaplan-meier",
		dist = "lnorm", add.line = TRUE, 
		include.cen = TRUE,  
		duplicate.points.method = "jitter",
		ylab = "log(manganese) log(ppb)", main = ""))
mtext(paste("Figure 15-4. Censored Probability Plot of",
	"Logged Manganese Sample", sep = "\n"), 
	cex = 1.25, font = 2, line = 1.75)

# Correlation between *UNIQUE* detected observations and censoring levels
#                     vs. plotting positions:
#------------------------------------------------------------------------
x <- qqplot.list$x
y <- qqplot.list$y

unique.y <- unique(y)
new.x <- x[match(unique.y, y)]

cor(new.x, unique.y)
#[1] 0.9887829



#-------------------------------------------
# Estimate mean and SD of log-transformed Mn
# based on Kaplan-Meier estimates
#-------------------------------------------

with(EPA.09.Ex.15.1.manganese.df,
	enparCensored(log(Manganese.ppb), Censored))

#Results of Distribution Parameter Estimation
#Based on Type I Censored Data
#--------------------------------------------
#
#Assumed Distribution:            None
#
#Censoring Side:                  left
#
#Censoring Level(s):              0.6931472 1.6094379 
#
#Estimated Parameter(s):          mean    = 2.3092890
#                                 sd      = 1.1816102
#                                 se.mean = 0.1682862
#
#Estimation Method:               Kaplan-Meier
#
#Data:                            log(Manganese.ppb)
#
#Censoring Variable:              Censored
#
#Sample Size:                     25
#
#Percent Censored:                24%

rm(qqplot.list, x, y, unique.y, new.x)


##################################################################################################


# Example 15-2, p. 15-16 to 15-20
#--------------------------------

# NOTE:  The Robust ROS method as presented in the EPA (2009) guidance 
#        is eqivalent to using the EnvStats functions enormCensored() or 
#        elnormCensored() (depending on whether you want to assume
#        a normal or lognormal distribution) and setting the arguments 
#        method="impute.w.qq.reg", prob.method="hirsch-stedinger", 
#        and plot.pos.con=0.  You can also set prob.method and plot.pos.con to 
#        other values and you would still be using the Robust ROS method.


# NOTE:  It is important to distinguish between the method used to compute plotting 
#        positions for a Q-Q plot versus the method used to estimate distribution 
#        parameters.  For the probability plots displayed in Figures 15-5 and 15-6, the 
#        method used to compute plotting positions is the Hirsh-Stedinger method with a 
#        plotting position constant of 0, and the method used to fit the line is 
#        least-squares.  So contrary to the titles that USEPA (2009) is using for these 
#        plots, these QQ-Plots are not based on the "Robust ROS" method.  The term 
#        "Robust ROS" refers to how the mean and standard deviation are estimated, 
#        not to how the Q-Q plot is created.

#####
# Raw Observations:

# Figure 15-5
#------------
windows()
qqplot.list <- with(EPA.09.Ex.15.1.manganese.df, 
	qqPlotCensored(Manganese.ppb, Censored, 
		prob.method = "hirsch-stedinger", plot.pos.con = 0,
		add.line = T,
		ylab = "Manganese Concentration (ppb)", main = ""))
mtext(paste("Figure 15-5. Robust ROS Censored Probability Plot",
    "of Manganese Concentrations", sep = "\n"), cex = 1.25, font = 2, 
    line = 1.75)

cbind(Detected.Value = qqplot.list$Ord[!qqplot.list$Censored], 
 Plotting.Position = round(qqplot.list$Cum[!qqplot.list$Censored], 3),
 Z.score = round(qqplot.list$x[!qqplot.list$Censored], 3))
#      Detected.Value Plotting.Position Z.score
# [1,]            3.3             0.245  -0.690
# [2,]            5.3             0.318  -0.474
# [3,]            6.3             0.356  -0.370
# [4,]            7.7             0.394  -0.270
# [5,]            8.4             0.432  -0.172
# [6,]            9.5             0.469  -0.077
# [7,]           10.0             0.507   0.018
# [8,]           11.9             0.545   0.114
# [9,]           12.1             0.583   0.210
#[10,]           12.6             0.621   0.308
#[11,]           16.9             0.659   0.410
#[12,]           17.9             0.697   0.515
#[13,]           21.6             0.735   0.627
#[14,]           22.7             0.773   0.748
#[15,]           34.5             0.811   0.880
#[16,]           45.9             0.848   1.030
#[17,]           53.6             0.886   1.207
#[18,]           77.2             0.924   1.434
#[19,]          106.3             0.962   1.776

# Correlation between *DETECTED* observations and plotting positions:
#-------------------------------------------------------------------
cor(qqplot.list$x[!qqplot.list$Censored], qqplot.list$y[!qqplot.list$Censored])
#[1] 0.9013241



#####
# Log-transformed Observations:

# Figure 15-6
#------------
windows()
qqplot.list <- with(EPA.09.Ex.15.1.manganese.df, 
	qqPlotCensored(Manganese.ppb, Censored, dist = "lnorm",
		prob.method = "hirsch-stedinger", plot.pos.con = 0,
		add.line = T, 
		ylab = "log(Manganese) log(ppb)", main = ""))
mtext(paste("Figure 15-6. Robust ROS Censored Probability Plot",
    "of Logged Manganese", sep = "\n"), line = 1.75, cex = 1.25, font = 2)

# Correlation between *DETECTED* observations and plotting positions:
#-------------------------------------------------------------------
cor(qqplot.list$x[!qqplot.list$Censored], qqplot.list$y[!qqplot.list$Censored])
#[1] 0.9943942

with(EPA.09.Ex.15.1.manganese.df, 
	elnormCensored(x = Manganese.ppb, censored = Censored, 
		method = "impute.w.qq.reg",
		prob.method = "hirsch-stedinger", plot.pos.con = 0))

#Results of Distribution Parameter Estimation
#Based on Type I Censored Data
#--------------------------------------------
#
#Assumed Distribution:            Lognormal
#
#Censoring Side:                  left
#
#Censoring Level(s):              2 5 
#
#Estimated Parameter(s):          meanlog = 2.277175
#                                 sdlog   = 1.261431
#
#Estimation Method:               impute.w.qq.reg
#
#Data:                            Manganese.ppb
#
#Censoring Variable:              Censored
#
#Sample Size:                     25
#
#Percent Censored:                24%


rm(qqplot.list)

##################################################################################################

# Example 15-3, p. 15-24
#-----------------------

Mn       <- EPA.09.Ex.15.1.manganese.df$Manganese.ppb
Censored <- EPA.09.Ex.15.1.manganese.df$Censored

Mn[Censored]
#[1] 5 2 5 5 2 2

# This example wants you to use a single censoring value.
# It appears that all observations <= 5 are set to censored. 
Censored <- Mn <= 5
Mn[Censored] <- 5

#Cohen's Method
#--------------
elnormCensored(x = Mn, censored = Censored, method = "mle")

#Results of Distribution Parameter Estimation
#Based on Type I Censored Data
#--------------------------------------------
#
#Assumed Distribution:            Lognormal
#
#Censoring Side:                  left
#
#Censoring Level(s):              5 
#
#Estimated Parameter(s):          meanlog = 2.323338
#                                 sdlog   = 1.201756
#
#Estimation Method:               mle
#
#Data:                            Mn
#
#Censoring Variable:              Censored
#
#Sample Size:                     25
#
#Percent Censored:                28%


#Parametric ROS
#--------------

elnormCensored(x = Mn, censored = Censored, 
	method = "qq.reg", plot.pos.con = 0.375)

#Results of Distribution Parameter Estimation
#Based on Type I Censored Data
#--------------------------------------------
#
#Assumed Distribution:            Lognormal
#
#Censoring Side:                  left
#
#Censoring Level(s):              5 
#
#Estimated Parameter(s):          meanlog = 2.325662
#                                 sdlog   = 1.250635
#
#Estimation Method:               qq.reg
#
#Data:                            Mn
#
#Censoring Variable:              Censored
#
#Sample Size:                     25
#
#Percent Censored:                28%


rm(Mn, Censored)

#####################################################################################################



