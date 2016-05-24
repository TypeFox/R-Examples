# Copyright (C) 2011 by 52 North Initiative for Geospatial Open Source Software GmbH, Contact: info@52north.org
# This program is free software; you can redistribute and/or modify it under the terms of the GNU General Public License version 2 as published by the Free Software Foundation. This program is distributed WITHOUT ANY WARRANTY; even without the implied WARRANTY OF MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program (see gpl-2.0.txt). If not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA or visit the Free Software Foundation web page, http://www.fsf.org.
# Author: Daniel Nuest (daniel.nuest@uni-muenster.de)
# Project: sos4R - visit the project web page, http://www.nordholmen.net/sos4r
library("sos4R");
library("ggplot2"); library("latticeExtra"); library("quantreg")

##############################################################################
# PegelOnlineSOS
pegelsos <- SOS(url = "http://pegelonline.wsv.de/webservices/gis/gdi-sos")

pegelsos

# what data do I get?
cat("\nNames of offerings:\n")
print(names(sosOfferings(pegelsos)))

# let's find interesting data
# Bake_Z: http://www.pegelonline.wsv.de/gast/stammdaten?pegelnr=9510066
procs <- sosProcedures(pegelsos)[["WASSERSTAND_ROHDATEN"]]
contain_bake <- procs %in% grep("*Wasserstand-Bake*", procs, value=TRUE)
baken <- subset(procs, contain_bake)
cat("\nGauges in Northern Sea:\n")
print(baken)

# what?
wasserstand_roh <- sosOfferings(pegelsos)[["WASSERSTAND_ROHDATEN"]]
# what?
wasserstand <- sosObservedProperties(wasserstand_roh)[1]
# PROBLEM: order of result list is not alswas the same!
# must be "Wasserstand":

# when?
tPeriod.days <- 3 # 30
tPeriod <- sosCreateEventTimeList(
		time = sosCreateTimePeriod(
				sos = pegelsos,
				begin = Sys.time() - (3600 * 24 * tPeriod.days),
				end = Sys.time()))
#tPeriod # str(tPeriod)
#encodeXML(tPeriod[[1]], pegelsos)

# three procedures, but only getting 1 element with one procedure...
cat("\nRequesting data... \n")
pegelObs <- getObservation(sos = pegelsos,
		observedProperty = wasserstand,
		offering = wasserstand_roh,
		procedure = baken[c(1, 3)],
		eventTime = tPeriod)

# show parts of the data frame:
cat("\nData excerpts:\n")
sosResult(pegelObs[[1]])[1:2,]
sosResult(pegelObs)[1:10,]

# not enough info? got field descriptions as attributes for each column:
cat("\nMetadata for two values:\n")
attributes(sosResult(pegelObs[[1]])$SamplingTime) # TIME
attributes(sosResult(pegelObs[[1]])$Wasserstand) # actual value

# do something with the data! (here we clean up first all values below 0)
r1 <- sosResult(pegelObs[[1]])
range(r1$Wasserstand)
r1clean <- subset(r1, Wasserstand > 0)
range(r1clean$Wasserstand)

r2 <- sosResult(pegelObs[[2]])
range(r2$Wasserstand)
r2clean <- subset(r2, Wasserstand > 0)
range(r2clean$Wasserstand)

plot(r1clean$SamplingTime, r1clean$Wasserstand, type = "l", ylim=c(200,800),
		lty = "solid", col = "blue",
		main = paste("Water level at", sosProcedures(pegelObs[c(1,2)])))
lines(r2clean$SamplingTime, r2clean$Wasserstand, type = "l", col = "orange",
		lty = "twodash",
		main = paste("Water level at", sosProcedures(pegelObs[[2]])))
par(ask = TRUE)

# Plot a quantile regression line with standard error bounds, using the quantreg package.
r1plot <- xyplot(r1clean$Wasserstand ~ r1clean$SamplingTime, r1clean, type = "l",
		col = "orange", main = paste0(sosProcedures(pegelObs[[1]]), 
				" with quantile regression line and error bounds"),
		xlab = "Time", ylab = "Water level")

r1plot <- r1plot + layer(panel.quantile(x, y, tau = c(.95, .5, .05)))
show(r1plot)

###################################
# Demo finished, try another one! #
###################################
