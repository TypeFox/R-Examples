#####
# File:     Chapter05.R
#
# Purpose:  Reproduce Examples in Chapter 5 of the book:
#
#           Millard, SP. (2013).  EnvStats: An R Package for 
#             Environmental Statistics.  Springer-Verlag.
#
# Author:   Steven P. Millard
#
# Last
# Updated:  10/14/2012
#####


#######################################################################################

library(EnvStats)

############################################################

# 5.2.1 Estimating Parameters of a Normal Distribution 
#-----------------------------------------------------

attach(EPA.94b.tccb.df)
enorm.list <- enorm(log(TcCB[Area == "Reference"]), ci = TRUE)
enorm.list


# Figure 5.1
#-----------
windows()
hist(log(TcCB[Area == "Reference"]), freq = FALSE, 
  xlim = c(-2, 1), xlab = "log [ TcCB (ppb) ]", 
  ylim = c(0, 1), col = "cyan", 
  main = "Log(TcCB) Data for Reference Area\nwith Fitted Normal Distribution")
params <- enorm.list$parameters
pdfPlot(dist = "norm", param.list = as.list(params), add = TRUE)

rm(enorm.list, params)
detach("EPA.94b.tccb.df")

###################################################################

# 5.2.2 Estimating Parameters of a Lognormal Distribution
#--------------------------------------------------------
attach(EPA.94b.tccb.df)
elnormAlt.list <- elnormAlt(TcCB[Area == "Reference"], ci = TRUE)
elnormAlt.list


# Figure 5.2
#-----------
windows()
hist(TcCB[Area == "Reference"], freq = FALSE, 
  xlim = c(0, 2), xlab = "TcCB (ppb)", 
  ylim = c(0, 2), col = "cyan", 
  main = "TcCB Data for Reference Area \nwith Fitted Lognormal Distribution")
params <- elnormAlt.list$parameters
pdfPlot(dist = "lnormAlt", param.list = as.list(params), add = TRUE)

rm(elnormAlt.list, params)
detach("EPA.94b.tccb.df")


sort(NIOSH.89.air.lead.vec)
round(elnormAlt(NIOSH.89.air.lead.vec, ci = TRUE)$interval$limits)
round(elnormAlt(NIOSH.89.air.lead.vec, ci = TRUE, ci.method = "zou")$interval$limits)

############################################################################

#5.3.1 Estimating Quantiles of a Normal Distribution
#----------------------------------------------------
EPA.09.Ex.21.1.aldicarb.df

attach(EPA.09.Ex.21.1.aldicarb.df)

eqnorm(Aldicarb.ppb[Well == "Well.1"], p = 0.95, ci = TRUE, ci.type = "lower", conf.level = 0.99)

sapply(split(Aldicarb.ppb, Well), function(x) {
  eqnorm(x, p = 0.95, ci = TRUE, ci.type = "lower", 
  conf.level = 0.99)$interval$limits["LCL"]
  })

sapply(split(Aldicarb.ppb, Well), function(x) {
  eqnorm(x, p = 0.95, ci = TRUE, ci.type = "upper", 
  conf.level = 0.99)$interval$limits["UCL"]
  })

detach("EPA.09.Ex.21.1.aldicarb.df")

####################################################################################

#5.3.2 Estimating Quantiles of a Lognormal Distribution
#------------------------------------------------------

EPA.09.Ex.17.3.chrysene.df

attach(EPA.09.Ex.17.3.chrysene.df)

Chrysene <- Chrysene.ppb[Well.type == "Background"]
eqlnorm(Chrysene, p = 0.95, ci = TRUE, ci.type = "upper", conf.level = 0.95)

rm(Chrysene)
detach("EPA.09.Ex.17.3.chrysene.df")

##################################################################################

#5.3.3 Estimating Quantiles of a Gamma Distribution
#--------------------------------------------------

attach(EPA.09.Ex.17.3.chrysene.df)

Chrysene <- Chrysene.ppb[Well.type == "Background"]

eqgamma(Chrysene, p = 0.95, ci = TRUE, ci.type = "upper", conf.level = 0.95)

rm(Chrysene)
detach("EPA.09.Ex.17.3.chrysene.df")

##################################################################################


#5.3.4 Nonparametric Estimates of Quantiles
#------------------------------------------

Sample.Size <- c(seq(5, 25, by = 5), 50, 75, 100)
conf.level <- tolIntNparConfLevel(Sample.Size, coverage = 0.95, ti.type = "upper")
cbind(Sample.Size, Confidence.Level = round(100 * conf.level))
rm(Sample.Size, conf.level)


EPA.09.Ex.21.6.nitrate.df[, 1:3]

EPA.09.Ex.21.6.nitrate.df[, c(2, 4:5)]

Nitrate <- EPA.09.Ex.21.6.nitrate.df$ Nitrate.mg.per.l
eqnpar(Nitrate, p = 0.95, ci = TRUE, ci.type = "lower", approx.conf.level = 0.95) 

eqnpar(Nitrate, p = 0.95, ci = TRUE, ci.type = "lower", lcl.rank = 10) 

eqnpar(Nitrate, p = 0.95, ci = TRUE, ci.type = "upper", approx.conf.level = 0.95) 

tolIntNparN(coverage = 0.95, ti.type = "upper", conf.level = 0.95)

rm(Nitrate)



