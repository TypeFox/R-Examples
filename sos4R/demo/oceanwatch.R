# Copyright (C) 2011 by 52 North Initiative for Geospatial Open Source Software GmbH, Contact: info@52north.org
# This program is free software; you can redistribute and/or modify it under the terms of the GNU General Public License version 2 as published by the Free Software Foundation. This program is distributed WITHOUT ANY WARRANTY; even without the implied WARRANTY OF MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program (see gpl-2.0.txt). If not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA or visit the Free Software Foundation web page, http://www.fsf.org.
# Author: Daniel Nuest (daniel.nuest@uni-muenster.de)
# Project: sos4R - visit the project web page, http://www.nordholmen.net/sos4r
library("sos4R")

##############################################################################
# Oceanwatch
# http://las.pfeg.noaa.gov/oceanWatch/oceanwatch.php
oceanwatch <- SOS(url = "http://oceanwatch.pfeg.noaa.gov/pysos/sos_mysql2.py",
		method = SosSupportedConnectionMethods()[["GET"]])
# warnings about missing response modes for offerings
warnings()

ocean.off <- sosOfferings(oceanwatch)
names(ocean.off)

ocean.proc <- sosProcedures(oceanwatch)
ocean.proc

ocean.obsProp <- sosObservedProperties(oceanwatch)
ocean.obsProp

sosCapabilitiesDocumentOriginal(oceanwatch)
# OK, but missing elements

describeSensor(oceanwatch, ocean.proc[[4]], inspect = TRUE)
# fails

lapply(X = ocean.off, FUN = getObservation, sos = oceanwatch)
# Es gab 36 Warnungen (Anzeige mit warnings())

length(ocean.off)
# 36

# -> SOS seems empty

# TODO

###################################
# Demo finished, try another one! #
###################################
