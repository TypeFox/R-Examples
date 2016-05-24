#####
# File:     Chapter01.R
#
# Purpose:  Reproduce Examples in Chapter 1 of the book:
#
#           Millard, SP. (2013).  EnvStats: An R Package for 
#             Environmental Statistics.  Springer-Verlag.
#
# Author:   Steven P. Millard
#
# Last
# Updated:  2012/12/13
#####

#######################################################################################

library(EnvStats)

############################################################


# 1.10.1 The TcCB data
#---------------------
EPA.94b.tccb.df

head(EPA.94b.tccb.df)

attach(EPA.94b.tccb.df)


############################################################


# 1.10.2 Computing summary statistics
#------------------------------------
summaryStats(TcCB ~ Area, data = EPA.94b.tccb.df)

summaryFull(TcCB ~ Area, data = EPA.94b.tccb.df)


#############################################################


# 1.10.3 Looking at the TcCB data
#--------------------------------

# Figure 1.1
#-----------
windows()
par(mar = c(4, 4, 4, 1) + 0.1)
stripChart(log(TcCB) ~ Area, data = EPA.94b.tccb.df, 
  col = c("red", "blue"), ylab = "Log  [ TcCB (ppb) ]")
mtext("Stripcharts of Log-Transformed TcCB Data", 
  line = 3, cex = 1.25, font = 2)


# Figure 1.2
#-----------
windows()
par(mfrow=c(2,1))
hist(log(TcCB[Area == "Reference"]), 
  xlim = c(-4, 6), col = "blue", 
  xlab = "log [ TcCB (ppb) ]", 
  ylab = "Number of Observations", 
  main = "Reference Area")
hist(log(TcCB[Area == "Cleanup"]), 
  xlim = c(-4, 6), nclass = 30, col = "red", 
  xlab = "log [ TcCB (ppb) ]", 
  ylab = "Number of Observations", 
  main = "Cleanup Area")


# Figure 1.3
#-----------
windows()
boxplot(log(TcCB) ~ Area, col = c("red", "blue"), 
  pars = list(outpch = 16), xlab = "Area", 
  ylab = "Log [ TcCB (ppb) ]",
  main = "Boxplots of Log-Transformed TcCB Data")


#############################################################


# 1.10.4 Quantile (Empirical CDF) Plots
#--------------------------------------

# Figure 1.4
#-----------
windows()
ecdfPlot(TcCB[Area=="Reference"], 
	xlab="TcCB (ppb)", 
	main = "Emipriacal CDF Plot of Reference Area TcCB")


# Figure 1.5
#-----------
windows()
cdfCompare(TcCB[Area=="Reference"], dist="lnorm", 
  xlab="Order Statistics for Reference Area TcCB (ppb)", 
  main=paste("Empirical CDF for ",
    "Reference Area TcCB (solid line)\n",
    "with Fitted Lognormal CDF (dashed line)", sep=""))
legend(0.65, 0.4, legend = c("Empirical CDF", "Fitted Lognormal CDF"), 
  lty = 1:2, col = c(4, 1), lwd = 3, bty = "n")


# Figure 1.6
#-----------
windows()
cdfCompare(log(TcCB[Area=="Reference"]), 
  log(TcCB[Area=="Cleanup"]), 
  xlab="Order Statistics for log [ TcCB (ppb) ]",
  main=paste("CDF for Reference Area (solid line)",
    "with CDF for Cleanup Area (dashed line)", 
    sep="\n"))
legend(1.5, 0.4, legend = c("Reference Area", "Cleanup Area"), lty = 1:2, 
  col = c(4, 1), lwd = 3, bty = "n")



##################################################################


# 1.10.5 Assessing goodness-of-fit with quantile-quantile plots
#--------------------------------------------------------------

# Figure 1.7
#-----------
windows()
qqPlot(TcCB[Area=="Reference"], dist="lnorm", 
  add.line=TRUE, points.col = "blue", 
  ylab="Quantiles of log [ TcCB (ppb) ]", 
  main=paste("Normal Q-Q Plot for",
    "Log-Transformed Reference Area TcCB Data", sep = "\n"))


# Figure 1.8
#-----------
windows()
qqPlot(TcCB[Area=="Reference"], dist="lnorm", plot.type="Tukey", 
  estimate.params=TRUE, add.line=TRUE, points.col = "blue", 
  main=paste("Tukey Mean-Difference Q-Q Plot for",
    "\nReference Area TcCB Fitted to ",
    "Lognormal Distribution", sep=""))


# Figure 1.9
#-----------
windows()
qqPlot(TcCB[Area=="Reference"], dist="gamma", 
  estimate.params = TRUE, digits = 2, 
  add.line=TRUE, points.col = "blue", 
  ylab="Quantiles of TcCB (ppb)", 
  main="Gamma Q-Q Plot for Reference Area TcCB Data")



####################################################################


# 1.10.6 Estimating distribution parameters
#------------------------------------------

elnorm(TcCB[Area=="Reference"], ci=TRUE)

elnormAlt(TcCB[Area=="Reference"], ci=TRUE)

egamma(TcCB[Area=="Reference"], ci=TRUE)

egammaAlt(TcCB[Area=="Reference"], ci=TRUE)



####################################################################


# 1.10.7 Testing for Goodness of Fit
#-----------------------------------

TcCB.ref <- TcCB[Area=="Reference"]

sw.lnorm <- gofTest(TcCB.ref, dist="lnorm")
sw.lnorm

# Figure 1.10
#------------
windows()
plot(sw.lnorm, digits=3)


sw.gamma <- gofTest(TcCB.ref, dist="gamma")
sw.gamma

# Figure 1.11
#------------
windows()
plot(sw.gamma, digits=3)


rm(TcCB.ref, sw.lnorm, sw.gamma)





####################################################################


# 1.10.8 Estimating Quantiles and Computing Confidence Limits
#------------------------------------------------------------

eqlnorm(TcCB[Area=="Reference"], p=0.9, ci=T)


####################################################################


# 1.10.9 Comparing Two Distributions Using Nonparametric Tests
#-------------------------------------------------------------
twoSampleLinearRankTest(TcCB[Area=="Cleanup"], 
	TcCB[Area=="Reference"], alternative="greater")

quantileTest(TcCB[Area=="Cleanup"], TcCB[Area=="Reference"], 
	alternative="greater", target.r=9)


####################################################################

detach("EPA.94b.tccb.df")

