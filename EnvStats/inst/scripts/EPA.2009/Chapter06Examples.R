#####
# File:     Script - EPA09, Chapter 06 Examples.R
#
# Purpose:  Reproduce Examples in Chapter 6 of the EPA Guidance Document
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
# Updated:  2013/10/15
#####

#NOTE:  Unless otherwise noted, differences between what is shown in the
#       EPA Guidance Document and results shown here are due to the 
#       EPA Guidance Document rounding intermediate results prior to 
#       the final result.

#######################################################################################

library(EnvStats)

# Figure 6-2, p. 6-15
#--------------------
windows()
plotPredIntNormTestPowerCurve(n = 10, k = 1, 
	range.delta.over.sigma = c(0, 5), 
	pi.type = "upper", conf.level = 0.99, 
	ylim = c(0, 1), 
	xlab = expression(paste(Delta, " (SDs above Background)")), 
	ylab = "Power", cex.main = 1.1, 
	main = "Figure 6-2. Normal Power Curve (n = 10) for 99% Prediction Limit Test")

##############################################################################

# Figure 6-3, p. 6-17
#--------------------
windows()
plotPredIntNormTestPowerCurve(n = 10, k = 1, 
	range.delta.over.sigma = c(0, 5), 
	pi.type = "upper", conf.level = 0.99, 
	ylim = c(0, 1), 
	xlab = "SD Units Above BG", 
	ylab = "Power",
	main = "Figure 6-3. EPA Reference Power Curves")

plotPredIntNormTestPowerCurve(n = 10, k = 2, conf.level = 0.99, 
	add = TRUE, plot.col = "red", plot.lty = 2)

plotPredIntNormTestPowerCurve(n = 10, k = 4, conf.level = 0.99, 
	add = TRUE, plot.col = "blue", plot.lty = 3)

legend("topleft", c("Quarterly", "Semi-Annual", "Annual"), lty = 3:1, 
	lwd = 3 * par("cex"), col = c("blue", "red", "black")) 

##################################################################

# Example 6-1, p. 6-18
#---------------------
predIntNormTestPower(n = 10, k = 1, 
	delta.over.sigma = (15 - 6) / 2, 
	pi.type = "upper", conf.level = 0.99)
#[1] 0.9012473



# Example 6-3, pp. 6-20 to 6-21
#------------------------------
stats   <- summaryStats(Sulfate.ppm ~ 1, data = EPA.09.Ex.6.3.sulfate.df,
	data.name = "Sulfate", digits = 1)
stats
#         N Mean   SD Median Min Max
#Sulfate 10  545 29.7    546 490 590

n     <- stats[, "N"]
SD    <- stats[, "SD"]
delta <- c(0, 200)

windows()
plotPredIntNormTestPowerCurve(n = n, k = 1, 
	axes = F,
	range.delta.over.sigma = range(delta / SD), 
	pi.type = "upper", conf.level = 0.99, 
	ylim = c(0, 1), 
	xlab = "Sulfate Conc. Increase (ppm)", 
	ylab = "Power",
	main = "Figure 6-3. Approximate s-Based Power Curve for Sulfate")
axis(2)
x.axis.ticks <- seq(0, 200, by = 50)
axis(1, at = x.axis.ticks / SD, labels = x.axis.ticks)
box()
abline(v = 75 / SD, lty = 2, lwd = 2)

predIntNormTestPower(n = n, k = 1, delta.over.sigma = 75 / SD, 
	pi.type = "upper", conf.level = 0.99)
#[1] 0.3916402

rm(stats, n, SD, delta, x.axis.ticks)
