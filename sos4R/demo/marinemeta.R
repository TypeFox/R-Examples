# Copyright (C) 2011 by 52 North Initiative for Geospatial Open Source Software GmbH, Contact: info@52north.org
# This program is free software; you can redistribute and/or modify it under the terms of the GNU General Public License version 2 as published by the Free Software Foundation. This program is distributed WITHOUT ANY WARRANTY; even without the implied WARRANTY OF MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program (see gpl-2.0.txt). If not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA or visit the Free Software Foundation web page, http://www.fsf.org.
# Author: Daniel Nuest (daniel.nuest@uni-muenster.de)
# Project: sos4R - visit the project web page, http://www.nordholmen.net/sos4r
library("sos4R")


##############################################################################
# OOSTethys SOS                                                              #
# http://www.oostethys.org/development/web-services/web-services-summary     #
##############################################################################
# Sensor Observation Service (SOS) for Marine Metadata Interoperability
# Initiative (MMI) # Using GET works!

# add converter for UOM C and others (which is not a valid unit)
myConverters <- SosDataFieldConvertingFunctions(
		# mapping for UOM:
		"C" = sosConvertDouble,
		"S/m" = sosConvertDouble,
		# mapping for definition:
		"http://mmisw.org/ont/cf/parameter/sea_water_salinity" = sosConvertDouble)
MBARI <- SOS("http://mmisw.org/oostethys/sos", method = "GET",
		dataFieldConverters = myConverters)

myOff <- sosOfferings(MBARI)[[1]]
myProc <- sosProcedures(MBARI)[[1]]

mbariObs <- getObservation(sos = MBARI, offering = myOff, procedure = myProc)

# explore data
data <- sosResult(mbariObs) # just one Observation in ObservationCollection, which is returned directly
cat("\nFirst lines of data:\n")
print(data[1:3,])
cat("\nMetadata attributes of one value:\n")
print(attributes(data$Temperature))
#str(data)

cat("\nCovariance of Temperature Conductivity Salinity:\n")
tcs.cov <- cov(data[, 5:7])
print(tcs.cov)
cat("\nCorrelation of Temperature Conductivity Salinity:\n")
tcs.cor <- cor(data[, 5:7])
print(tcs.cor)

# quick plot
plot(data)
par(ask = TRUE)
plot(data[,c("Temperature", "Salinity", "Conductivity")],
		main = "MBARI MMI SOS",
		sub = "OOSTethys Sensor Observation Service")

cat("\nExpect a warning here because of incomplete swe:Quantitiy.\n")

###################################
# Demo finished, try another one! #
###################################
