#####
# File:     Chapter08.R
#
# Purpose:  Reproduce Examples in Chapter 8 of the book:
#
#           Millard, SP. (2013).  EnvStats: An R Package for 
#             Environmental Statistics.  Springer-Verlag.
#
# Author:   Steven P. Millard
#
# Last
# Updated:  2013/08/27
#####


#######################################################################################

library(EnvStats)

############################################################

# 8.4.1 Quantile (Empirical CDF) Plots for Censored Data 
#-----------------------------------------------------

# Figures 8.1 & 8.2
#------------------

windows()
with(Helsel.Cohn.88.silver.df,
	ecdfPlotCensored(Ag, Censored, 
	xlab = expression(paste("Ag (", mu, "g/L)", sep = "")),
	main = "Empirical CDF Plot for Silver\n(Method of Michael-Schucany)", 
	include.cen = TRUE))

windows()
with(Helsel.Cohn.88.silver.df,
	ecdfPlotCensored(log.Ag, Censored, 
	xlab = expression(paste("log [ Ag (", mu, "g/L) ]", sep = "")), 
	main = "Empirical CDF Plot for Silver\n(Method of Michael-Schucany)", 
	include.cen = TRUE))


#######################################################################################


# 8.4.2 Comparing an Empirical and Hypothesized CDF
#--------------------------------------------------

# Figures 8.3 & 8.4
#------------------

windows()
with(Helsel.Cohn.88.silver.df,
	cdfCompareCensored(Ag, Censored, distribution = "lnorm", 
	xlab = expression(paste("Ag (", mu, "g/L)", sep = ""))))

windows()
with(Helsel.Cohn.88.silver.df,
	cdfCompareCensored(log(Ag), Censored, distribution="norm", 
	xlab = expression(paste("log [ Ag (", mu, "g/L) ]", sep = ""))))


#########################################################################################

# 8.4.3 Comparing Two Empirical CDFs
#-----------------------------------

# Figure 8.5
#-----------

attach(Millard.Deverel.88.df)
Cu.Alluvial <- Cu[Zone == "Alluvial.Fan"]
Cu.Alluvial.cen <- Cu.censored[Zone == "Alluvial.Fan"]
Cu.Basin <- Cu[Zone == "Basin.Trough"]
Cu.Basin.cen <- Cu.censored[Zone == "Basin.Trough"]

windows()
cdfCompareCensored(x = Cu.Alluvial, censored = Cu.Alluvial.cen, 
	y = Cu.Basin, y.censored = Cu.Basin.cen)

detach("Millard.Deverel.88.df")
rm(Cu.Alluvial, Cu.Alluvial.cen, Cu.Basin, Cu.Basin.cen)

#####################################################################################

# 8.4.4 Q-Q Plots for Censored Data
#----------------------------------

# Figure 8.6
#-----------

windows()
with(Modified.TcCB.df, 
	qqPlotCensored(TcCB, Censored, distribution = "lnorm", 
	add.line = TRUE, points.col = "blue"))

###################################################################################

# 8.4.5 Box-Cox Transformations for Censored Data
#------------------------------------------------

# Figure 8.7
#-----------

boxcox.list <- with(Modified.TcCB.df, 
	boxcoxCensored(TcCB, Censored))

windows()
plot(boxcox.list)

rm(boxcox.list)

#################################################################################

# 8.5 Estimating Distribution Parameters
#---------------------------------------

# 8.5.1 The Normal Distribution
#------------------------------

# Estimating Parameters of a Normal Distribution
#-----------------------------------------------

with(EPA.09.Ex.15.1.manganese.df, 
	elnormCensored(Manganese.ppb, Censored))

with(EPA.09.Ex.15.1.manganese.df, 
	elnormCensored(Manganese.ppb, Censored, method = "qq.reg"))

with(EPA.09.Ex.15.1.manganese.df, 
	elnormCensored(Manganese.ppb, Censored, 
	method = "impute.w.qq.reg", prob.method = "hirsch-stedinger"))

with(EPA.09.Ex.15.1.manganese.df, 
	elnormCensored(Manganese.ppb, Censored, 
	method = "impute.w.qq.reg", prob.method = "hirsch-stedinger",
	plot.pos.con = 0))



# Confidence Intervals for the Mean of the Normal Distribution
#-------------------------------------------------------------

with(EPA.09.Ex.15.1.manganese.df, 
	elnormCensored(Manganese.ppb, Censored, ci = TRUE))


#################################################################################

# 8.5.2 The Lognormal Distribution
#---------------------------------

with(Modified.TcCB.df, 
	elnormAltCensored(TcCB, Censored, method = "mle", ci = TRUE))

#################################################################################

# 8.5.3 The Gamma Distribution
#-----------------------------

with(Modified.TcCB.df, 
	egammaAltCensored(TcCB, Censored, ci = TRUE, ci.method = "profile.likelihood"))


#################################################################################

# 8.5.4 Estimating the Mean Nonparametrically
#--------------------------------------------

with(EPA.09.Ex.15.1.manganese.df,
	enparCensored(log(Manganese.ppb), Censored, ci = TRUE))

with(EPA.09.Ex.15.1.manganese.df,
	enparCensored(Manganese.ppb, Censored, ci = TRUE))


#################################################################################

# 8.6 Estimating Distribution Quantiles
#--------------------------------------

# 8.6.1 Parametric Estimates of Quantiles
#----------------------------------------

#--------------------------------------
# The Normal and Lognormal Distribution
#--------------------------------------

# First compute the estimated 90th percentile and 
# associated 95% confidence interval using the complete data.
#------------------------------------------------------------
with(EPA.94b.tccb.df, eqlnorm(TcCB[Area == "Reference"], p = 0.9, ci = TRUE))

# Now compute the estimated 90th percentile and confidence interval 
# using the censored data
#-------------------------------------------------------------------

# Estimate the mean and stanrdard deviation using the MLE method,
# and compute the confidence interval using the formula as if 
# we had complete data.
#-----------------------------------------------------------------
with(Modified.TcCB.df, eqlnormCensored(TcCB, Censored, p = 0.9, ci = TRUE))


# Estimate the mean and stanrdard deviation using imputation 
# with Q-Q regression (robust ROS), and compute the confidence 
# interval using the gpq method.
#-----------------------------------------------------------------

with(Modified.TcCB.df, eqlnormCensored(TcCB, Censored, p = 0.9, 
	method = "impute.w.qq.reg", ci = TRUE, 
	ci.method = "gpq", seed = 47))



#--------------------
# Other Distributions
#--------------------

egamma.list <- with(Modified.TcCB.df, egammaCensored(TcCB, Censored))

eqgamma(egamma.list, p = 0.9)


# Clean up
#---------
rm(egamma.list)

 
##########################################################################

# 8.7 Prediction Intervals
#-------------------------

# 8.7.1 Parametric Prediction Intervals
#--------------------------------------

enorm.list <- with(EPA.09.Ex.18.3.TCE.df, 
	enormCensored(TCE.ppb[Well.type == "Background"], 
		Censored[Well.type == "Background"]))

predIntNorm(enorm.list, k = 4, pi.type = "upper", conf.level = 0.95, method = "exact")

# Clean up
#---------
rm(enorm.list)

##############################################################################

# 8.8 Tolerance Intervals
#------------------------

# 8.8.1 Parametric Tolerance Intervals
#-------------------------------------

# Estimate the mean and stanrdard deviation using the MLE method,
# and compute the tolerance interval using the formula as if 
# we had complete data.
#-----------------------------------------------------------------
with(Modified.TcCB.df, tolIntLnormCensored(TcCB, Censored, ti.type = "upper"))


# Estimate the mean and stanrdard deviation using imputation 
# with Q-Q regression (robust ROS), and compute the tolerence 
# interval using the gpq method.
#-----------------------------------------------------------------
with(Modified.TcCB.df, tolIntLnormCensored(TcCB, Censored, ti.type = "upper", 
	method = "impute.w.qq.reg", ti.method = "gpq", seed = 47))


#############################################################################

# 8.9 Hypothesis Tests
#---------------------

# 8.9.1 Goodness-of-Fit Tests
#----------------------------

sw.list.norm <- with(Modified.TcCB.df, gofTestCensored(TcCB, Censored))
sw.list.norm

sw.list.lnorm <- with(Modified.TcCB.df, gofTestCensored(TcCB, Censored, dist = "lnorm"))
sw.list.lnorm

windows()
plot(sw.list.norm)

windows()
plot(sw.list.lnorm)

# Clean up
#---------
rm(sw.list.norm, sw.list.lnorm)

##############################################################################

# 8.9.2 Nonparametric Tests to Compare Two Groups
#------------------------------------------------

attach(Millard.Deverel.88.df)

Cu.AF <- Cu[Zone=="Alluvial.Fan"]
Cu.AF.cen <- Cu.censored[Zone=="Alluvial.Fan"]
Cu.BT <- Cu[Zone=="Basin.Trough"]
Cu.BT.cen <- Cu.censored[Zone=="Basin.Trough"]
twoSampleLinearRankTestCensored(Cu.AF, Cu.AF.cen, Cu.BT, Cu.BT.cen, 
	test = "normal.scores.2", var = "hypergeometric")

# Clean up
#---------
rm(Cu.AF, Cu.AF.cen, Cu.BT, Cu.BT.cen)

detach("Millard.Deverel.88.df")



