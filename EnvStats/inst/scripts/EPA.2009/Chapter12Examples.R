#####
# File:     Script - EPA09, Chapter 12 Examples.R
#
# Purpose:  Reproduce Examples in Chapter 12 of the EPA Guidance Document
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
# Updated:  January 15, 2012
#####

#NOTE:  Unless otherwise noted, differences between what is shown in the
#       EPA Guidance Document and results shown here are due to the 
#       EPA Guidance Document rounding intermediate results prior to 
#       the final result.

#######################################################################################

library(EnvStats)

# Example 12-1, pp. 12-3 to 12-5
#-------------------------------

#NOTE:  It is not clear what plotting position constant is 
#       being used.  In EnvironmentalStats for R, the default
#       value of the plotting position is 0.375 for an 
#       assumed normal or lognormal distribution.
#
#       The EPA Errata Sheet gives updated values for correlations for the
#       probability plots.

CCL4 <- EPA.09.Ex.12.1.ccl4.df$CCL4.ppb

#####
# Figure 12-1

windows()
dum.list <- qqPlot(CCL4, add.line = T,
	ylim = c(0, 8000), ylab = "Carbon tetrachloride (ppb)", 
	main = "")
Cor <- round(cor(dum.list$x, dum.list$y), 3)
title(main = bquote(paste("Figure 12-1. Probability Plot on Raw Concentrations (",
	italic(r) == .(Cor), ")", sep = "")))


#####
# Figure 12-2

windows()
dum.list <- qqPlot(log(CCL4), add.line = T,
	ylim = c(0, 10), ylab = "Log(carbon tetrachloride) log(ppb)", 
	main = "")
Cor <- round(cor(dum.list$x, dum.list$y), 3)
title(main = bquote(paste("Figure 12-2. Probability Plot on Logged Observations (",
	italic(r) == .(Cor), ")", sep = "")))

# Or alternatively:

dum.list <- qqPlot(CCL4, dist = "lnorm", add.line = T, 
	ylim = c(0, 10), ylab = "Log(carbon tetrachloride) log(ppb)", 
	main = "")
Cor <- round(cor(dum.list$x, dum.list$y), 3)
title(main = bquote(paste("Figure 12-2. Probability Plot on Logged Observations (",
	italic(r) == .(Cor), ")", sep = "")))




#####
# Remove "outlier"

CCL4.w.out.Max <- CCL4[CCL4 != max(CCL4)]


#####
# Figure 12-3

windows()
dum.list <- qqPlot(CCL4.w.out.Max, add.line = T,
	ylim = c(0, 400), ylab = "Carbon tetrachloride (ppb)", 
	main = "")
Cor <- round(cor(dum.list$x, dum.list$y), 3)
title(main = bquote(paste("Figure 12-3. Outlier-Deleted Probability Plot on Original Scale (",
	italic(r) == .(Cor), ")", sep = "")))



#####
# Figure 12-4

windows()
dum.list <- qqPlot(CCL4.w.out.Max, dist = "lnorm", add.line = T, 
	ylim = c(0, 10), ylab = "Log(carbon tetrachloride) log(ppb)", 
	main = "")
Cor <- round(cor(dum.list$x, dum.list$y), 3)
title(main = bquote(paste("Figure 12-4. Outlier-Deleted Probability Plot on Logarithmic Scale (",
	italic(r) == .(Cor), ")", sep = "")))



rm(CCL4, CCL4.w.out.Max, dum.list, Cor)


#######################################################################################


# Example 12-2, pp. 12-6 to 12-7
#-------------------------------

CCL4 <- EPA.09.Ex.12.1.ccl4.df$CCL4.ppb

windows()
par(mfrow = c(1,2), mar = c(1, 4, 1, 1), oma = c(0, 0, 3, 0))
boxplot(CCL4, log = "y", ylab = "carbon tetrachloride (ppb)")
points(1, mean(CCL4), pch = 3)
boxplot(log(CCL4), ylab = "log(carbon tetrachloride) log(ppb)")
points(1, mean(log(CCL4)), pch = 3)
mtext("Figure 12-5. Comparative Carbon Tetrachloride Boxplots Indicating Outliers", 
	line = 1, cex = 1.25, outer = T)

rm(CCL4)


########################################################################################


# Example 12-3, pp. 12-9 to 12-10
#--------------------------------

#NOTE:  you need to install the package "outliers" in order to 
#       perform Dixon's test for outliers

library(outliers)

CCL4 <- EPA.09.Ex.12.1.ccl4.df$CCL4.ppb

dixon.test(log(CCL4), type = 0, opposite = FALSE, two.sided = FALSE)

#Results of Hypothesis Test
#--------------------------
#
#Alternative Hypothesis:          highest value 8.8630498279191 is an outlier
#
#Test Name:                       Dixon test for outliers
#
#Data:                            log(CCL4)
#
#Test Statistic:                  Q = 0.4509385
#
#P-value:                         0.04928552

rm(CCL4)


#####################################################################################


# Example 12-4, pp. 12-12 to 12-14
#---------------------------------


#NOTE:  It is not clear what plotting position constant is 
#       being used.  In EnvironmentalStats for R, the default
#       value of the plotting position is 0.375 for an 
#       assumed normal or lognormal distribution.

Naphthalene <- EPA.09.Ex.12.4.naphthalene.df$Naphthalene.ppb

#####
# Figure 12-6

windows()
dum.list <- qqPlot(Naphthalene, add.line = T,
	ylim = c(0, 40), 
	ylab = "Naphthalene (ppb)", 
	main = "Figure 12-6. Naphthalene Probability Plot")
round(cor(dum.list$x, dum.list$y), 3)
#[1] 0.735


#####
# Figure 12-7

windows()
dum.list <- qqPlot(Naphthalene, dist = "lnorm", add.line = T, 
	ylim = c(0, 4), ylab = "Log Naphthalene log(ppb)", 
	main = "Figure 12-7. Log Naphthalene Probability Plot")
round(cor(dum.list$x, dum.list$y), 3)
#[1] 0.949

#NOTE:  Rosner's test for outliers is not available in R

rm(Naphthalene, dum.list)

