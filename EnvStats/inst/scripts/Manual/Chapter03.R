#####
# File:     Chapter03.R
#
# Purpose:  Reproduce Examples in Chapter 3 of the book:
#
#           Millard, SP. (2013).  EnvStats: An R Package for 
#             Environmental Statistics.  Springer-Verlag.
#
# Author:   Steven P. Millard
#
# Last
# Updated:  12/14/2012
#####


#######################################################################################

library(EnvStats)

############################################################

# 3.3.1 Summary Statistics for TcCB Concentrations
#-------------------------------------------------
EPA.94b.tccb.df


# 3.4 Empirical PDF Plots
#------------------------

# Figure 3.1
#-----------
attach(EPA.94b.tccb.df)

log.TcCB <- log(TcCB[Area == "Reference"])

windows()
hist(log.TcCB, freq = FALSE,  xlim = c(-2, 1),  
  col = "grey",  xlab = "log [ TcCB (ppb) ]", 
  ylab = "Relative Frequency", 
  main = "Histogram and Empirical PDF Plot for\n Log-Transformed Reference Area TcCB Data")
epdfPlot(log.TcCB, epdf.col = "blue", add = TRUE)

rm(log.TcCB)
detach("EPA.94b.tccb.df")


# 3.6.1 Q-Q Plots for the Normal and Lognormal Distribution
#----------------------------------------------------------

# Figure 3.2
#-----------
attach(EPA.94b.tccb.df)

windows()
qqPlot(TcCB[Area == "Reference"], add.line = TRUE, 
  points.col = "blue", ylab = "Quantiles of TcCB (ppb)", 
  main = "Normal Q-Q Plot for Reference Area TcCB Data")

detach("EPA.94b.tccb.df")



# 3.6.2 Q-Q Plots for Other Distributions
#----------------------------------------
EPA.92c.benzene1.df

# Figure 3.3
#-----------
attach(EPA.92c.benzene1.df)

Benzene[Censored] <- 1

windows()
qqPlot(Benzene, dist = "pois", estimate.params = TRUE, 
  duplicate.points.method = "number", add.line=TRUE, 
  qq.line.type = "0-1", points.col = "blue")

# Figure 3.4
#-----------
windows()
set.seed(721)
qqPlot(Benzene, dist = "pois", estimate.params = TRUE, 
  duplicate.points.method = "jitter", add.line = TRUE, 
  qq.line.type = "0-1", points.col = "blue")

rm(Benzene)
detach("EPA.92c.benzene1.df")


# 3.6.3 Using Q-Q Plots to Compare Two Data Sets
#-----------------------------------------------

# Figure 3.5
#-----------
attach(EPA.94b.tccb.df)
windows()
qqPlot(log(TcCB[Area == "Reference"]), log(TcCB[Area == "Cleanup"]), 
  points.col = "blue", plot.pos.con = 0.375, equal.axes = TRUE, 
  add.line = TRUE, qq.line.type="0-1", 
  xlab = paste("Quantiles of log [ TcCB (ppb) ]",
    "for Reference Area"), 
  ylab = paste("Quantiles of log [ TcCB (ppb) ]",
    "for Cleanup Area"), 
  main = paste("Q-Q Plot Comparing Cleanup and",
    "\nReference Area TcCB Data"))
detach('EPA.94b.tccb.df')


# 3.6.4 Building an Internal Gestalt for Q-Q Plots
#-------------------------------------------------

# Figure 3.6
#-----------
set.seed(426)
windows()
qqPlotGestalt(num.pages = 1, add.line = TRUE, ask = FALSE, 
  points.col = "blue")


# Figure 3.7
#-----------
windows()
qqPlotGestalt(num.pages = 1, add.line = TRUE, ask = FALSE, 
  plot.type = "Tukey", estimate.params = TRUE, 
  points.col = "blue")



# 3.7 Box-Cox Data Transformations and Q-Q Plots
#-----------------------------------------------

# Figure 3.8
#-----------
attach(EPA.94b.tccb.df)
boxcox.list <- boxcox(TcCB[Area == "Reference"])
windows()
plot(boxcox.list)

plot(boxcox.list, plot.type = "Q-Q")

plot(boxcox.list, plot.type = "Tukey")

plot(boxcox.list, plot.type = "All")

boxcox.list

rm(boxcox.list)
detach('EPA.94b.tccb.df')

