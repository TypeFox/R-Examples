# Copyright (C) 2010 by 52 North Initiative for Geospatial Open Source Software GmbH, Contact: info@52north.org
# This program is free software; you can redistribute and/or modify it under the terms of the GNU General Public License version 2 as published by the Free Software Foundation. This program is distributed WITHOUT ANY WARRANTY; even without the implied WARRANTY OF MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program (see gpl-2.0.txt). If not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA or visit the Free Software Foundation web page, http://www.fsf.org.
# Author: Daniel Nuest (daniel.nuest@uni-muenster.de)
# Project: sos4R - visit the project web page, http://www.nordholmen.net/sos4r
library("sos4R")
sessionInfo()

# establish a connection to a SOS instance with default settings
weathersos <- SOS(url = 
				"http://v-swe.uni-muenster.de:8080/WeatherSOS/sos")
weathersos
summary(weathersos)

################################################################################
# Request data for specific time intervall and offering
names(sosOfferings(weathersos))
sosOfferings(weathersos)[[1]]
off <- sosOfferings(weathersos)[["ATMOSPHERIC_TEMPERATURE"]]
sosObservedProperties(off)
# Exploration of new services is NOT the strong side of sos4R

# offering is only mandatory parameter, but use time interval to limit data
obs <- getObservation(sos = weathersos, offering = off,
							# alternative: offering = "ATMOSPHERIC_TEMPERATURE"
		eventTime = sosCreateTime(sos = weathersos,
				time = "2009-08-10 12:00/2009-08-20 12:00"))
obs
#length(obs)
#str(obs)

# subsettable in the normal way
sosProcedures(obs[1])
sosProcedures(obs[[2]])

# accessor functions
sosBoundedBy(obs)
sosCoordinates(obs[[2]])

sosProcedures(obs)
# vs.
sosProcedures(weathersos)

################################################################################
# get the REAL data
result <- sosResult(obs) # sosResult(obs[[1]])
summary(result)
dim(result); dim(sosResult(obs[[2]]))

# shortcut to spatial data
obs.spdf <- as(obs, "SpatialPointsDataFrame")
summary(obs.spdf)
# Metadata are useful!

################################################################################
# Create plot
# Attention: plots ignore the fact that the times do NOT perfectly match!
x <- 800
plot(x = obs[[1]]@result[[1]][1:x], y = obs[[1]]@result[[3]][1:x],
		type = "l", lwd = "2",
		col = "steelblue", main = "Temperature in Muenster and Kaernten, 2009",
		xlab = "Time (00:00 o'clock)",
		ylab = "Temperature (degree C)",
		xaxt="n") # do not plot x-axis
r <- as.POSIXct(round(range(obs[[1]]@result[[1]]), "days"))
axis.POSIXct(side = 1, x = obs[[1]]@result[[1]][1:x], format = "%d. %h",
		at = seq(r[1], r[2], by="days"))
lines(x = obs[[2]]@result[[1]][1:x], y = obs[[2]]@result[[3]][1:x],
		col = "orange", lwd = "2")
legend("topleft", legend = c("Muenster", "Kaernten"),
		col = c("steelblue", "orange"), lty = 1, bty="n")

################################################################################
# Time series analysis example
timeSeriesDemo()

# Is sos4R useful to you?
timeSeriesDemo(inspect = TRUE)

timeSeriesDemo <- function(inspect = FALSE) {
	library("xts")
	
	# set up request parameters
	station <- sosProcedures(weathersos)[[1]]
	temperatureOffering <- sosOfferings(weathersos)[["ATMOSPHERIC_TEMPERATURE"]]
	temperature <- sosObservedProperties(temperatureOffering)[1]
	lastWeek <- sosCreateTimePeriod(sos = weathersos,
			begin = as.POSIXct(Sys.time() - 3600 * 24 * 7),
			end = as.POSIXct(Sys.time()))
	
	# make the request
	obsLastWeek <- getObservation(sos = weathersos, # verbose = TRUE,
			observedProperty = temperature,
			procedure = station,
			eventTime = sosCreateEventTimeList(lastWeek),
			offering = temperatureOffering,
			inspect = inspect)
	dataLastWeek <- sosResult(obsLastWeek)
	
	# inspect data
	#summary(dataLastWeek)
	#names(dataLastWeek)
	
	# create time series from data and plot
	tempLastWeek <- xts(x = dataLastWeek[["urn:ogc:def:property:OGC::Temperature"]],
			order.by = dataLastWeek[["Time"]])
	
	# calculate regression (polynomial fitting)
	temp <- dataLastWeek[["urn:ogc:def:property:OGC::Temperature"]]
	time <- as.numeric(dataLastWeek[["Time"]])
	x = loess(temp~time,
			na.omit(dataLastWeek),enp.target=10)
	
	# create plot
	plot(tempLastWeek, main = "Temperature at WeatherSOS-Station in Muenster - LAST WEEK",
			xlab = "Time",
			ylab = paste("Temperature in", attributes(temp)[["unit of measurement"]]),
			major.ticks = "days")
	lines(dataLastWeek$Time, x$fitted, col = 'red', lwd=3)
	
	cat("data values: ")
	print(periodicity(tempLastWeek))
}

################################################################################
# Spatial analysis example, running while taking questions
demo("airquality")

