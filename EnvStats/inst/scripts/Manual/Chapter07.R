#####
# File:     Chapter07.R
#
# Purpose:  Reproduce Examples in Chapter 7 of the book:
#
#           Millard, SP. (2013).  EnvStats: An R Package for 
#             Environmental Statistics.  Springer-Verlag.
#
# Author:   Steven P. Millard
#
# Last
# Updated:  2013/06/09
#####


#######################################################################################

library(EnvStats)

############################################################

# 7.2.1 One Sample Goodness-of-Fit Tests for Normality 
#-----------------------------------------------------

attach(EPA.94b.tccb.df)
TcCB.Ref <- TcCB[Area == "Reference"]
sw.list.norm <- gofTest(TcCB.Ref)
sw.list.norm

sw.list.lnormAlt <- gofTest(TcCB.Ref, dist = "lnormAlt")
sw.list.lnormAlt

# Figure 7.1
#-----------
windows()
plot(sw.list.norm)

# Figure 7.2
#-----------
windows()
plot(sw.list.lnormAlt)

windows()
plot(sw.list.norm, plot.type = "Tukey")

detach("EPA.94b.tccb.df")
rm(TcCB.Ref, sw.list.norm, sw.list.lnormAlt)

#######################################################################################

# 7.2.2 Testing Several Groups for Normality
#-------------------------------------------

EPA.09.Ex.10.1.nickel.df

# Figure 7.3
#-----------
windows()
stripChart(Nickel.ppb ~ Well, data = EPA.09.Ex.10.1.nickel.df, 
	show.ci = FALSE, points.cex = 1.5, col = 1:4, ylab = "Nickel (ppb)")
mtext("Nickel Concentrations by Well", line = 3, font = 2, cex = 1.5)

# Figure 7.4
#-----------
windows()
stripChart(log10(Nickel.ppb) ~ Well, data = EPA.09.Ex.10.1.nickel.df, 
	show.ci = FALSE, points.cex = 1.5, col = 1:4, 
	ylab = expression(paste(log[10], "[ Nickel (ppb) ]")))
mtext("Log10-Transformed Nickel Concentrations by Well", line = 3, font = 2, cex = 1.5)


sw.list <- gofGroupTest(Nickel.ppb ~ Well, 
	data = EPA.09.Ex.10.1.nickel.df)
sw.list

# Figure 7.5
#-----------
windows()
plot(sw.list)


sw.log.list <- gofGroupTest(Nickel.ppb ~ Well, 
	data = EPA.09.Ex.10.1.nickel.df, dist = "lnorm")
sw.log.list

# Figure 7.6
#-----------
windows()
plot(sw.log.list)

rm(sw.list, sw.log.list)


###################################################################################################

# 7.2.3 One Sample Goodness-of-Fit Tests for Other Distributions
#---------------------------------------------------------------

attach(EPA.92c.benzene1.df)

Benzene[Censored] <- 1

table(Benzene)

lambda.hat <- epois(Benzene)$parameters

lambda.hat

ppois(c(-1, 0, 2, Inf), lambda = lambda.hat)

diff(ppois(c(-1, 0, 2, Inf), lambda = lambda.hat))

chisq.list <- gofTest(Benzene, test = "chisq", dist = "pois", cut.points = c(-1, 0, 2, Inf))

# Figure 7.7
#-----------
windows()
plot(chisq.list)


rm(Benzene, lambda.hat, chisq.list)
detach("EPA.92c.benzene1.df")


########################################################################################################

# 7.2.4 Two-Sample Goodness-of-Fit Test to Compare Samples
#---------------------------------------------------------

attach(EPA.94b.tccb.df)

log.TcCB.ref <- log(TcCB[Area=="Reference"])

log.TcCB.clean <- log(TcCB[Area=="Cleanup"])

ks.list <- gofTest(x = log.TcCB.ref, y = log.TcCB.clean, test = "ks")

ks.list

# Figure 7.8
#-----------
windows()
plot(ks.list)

rm(log.TcCB.ref, log.TcCB.clean, ks.list)
detach("EPA.94b.tccb.df")

#########################################################################################################

# 7.3.1 Two- and k-Sample Comparisons for Location
#-------------------------------------------------
summaryStats(log(TcCB) ~ Area, data = EPA.94b.tccb.df, 
	digits = 2, p.value = TRUE, stats.in.rows = TRUE)


summaryStats(log(TcCB) ~ Area, data = EPA.94b.tccb.df, 
	digits = 2, p.value = TRUE, stats.in.rows = TRUE,
	test.arg.list = list(var.equal = FALSE))

windows()
stripChart(log(TcCB) ~ Area, data = EPA.94b.tccb.df, 
	col = c("red", "blue"), p.value = TRUE,
	ylab = "Log [ TcCB (ppb) ]")
mtext("Strip Chart of TcCB Data with Test Results", cex = 1.5, line = 3, font = 2)

#################################################################################################

# 7.3.2 Chen’s Modified One-Sample t-Test for Skewed Data
#--------------------------------------------------------
hist(EPA.02d.Ex.9.mg.per.L.vec, col = "cyan", xlab = "Concentration (mg/L)", 
	main = "Concentrations at Exposure Unit for Exhibit 9 (USEPA, 2002d)")

gofTest(EPA.02d.Ex.9.mg.per.L.vec)$p.value

gofTest(EPA.02d.Ex.9.mg.per.L.vec, dist = "lnorm")$p.value

chenTTest(EPA.02d.Ex.9.mg.per.L.vec, mu = 50)


#################################################################################################

# 7.5.1 Testing for Trend in the Presence of Seasons
#---------------------------------------------------
Total.P.df

with(Total.P.df, plot(CB3.3e, type = "o", xaxt = "n", 
	xlab = "Time", ylab = "Total P (mg)"))

with(Total.P.df, axis(1, at = (1:length(CB3.3e))[Month == "Jan"], 
	labels = Year[Month == "Jan"]))

title(main = "Total Phosphorus at Station CB3.3e")

kendallSeasonalTrendTest(CB3.3e ~ Month + Year, data = Total.P.df)






