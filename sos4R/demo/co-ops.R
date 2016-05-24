# Copyright (C) 2011 by 52 North Initiative for Geospatial Open Source Software GmbH, Contact: info@52north.org
# This program is free software; you can redistribute and/or modify it under the terms of the GNU General Public License version 2 as published by the Free Software Foundation. This program is distributed WITHOUT ANY WARRANTY; even without the implied WARRANTY OF MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program (see gpl-2.0.txt). If not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA or visit the Free Software Foundation web page, http://www.fsf.org.
# Author: Daniel Nuest (daniel.nuest@uni-muenster.de)
# Project: sos4R - visit the project web page, http://www.nordholmen.net/sos4r
library("sos4R")


################################################################################
# NOAA Center for Operational Oceanographic Products and Services (CO-OPS)
# http://oceanservice.noaa.gov/programs/coops/
#
# CO-OPS SOS Services
#
# Related to demo "IOOS", CO-OPS' Implementation of IOOS
#
# https://geossregistries.info/geosspub/component_details_ns.jsp?compId=urn:uuid:c1af67f9-4a1b-42d2-b352-e2fdb3bcdeb1
# request examples at http://opendap.co-ops.nos.noaa.gov/ioos-dif-sos/
# Also good information about existing/active stations!
#
# Presentation: http://sdf.ndbc.noaa.gov/sos/IOOS_DIF_SOS_Project.ppt
#

################################################################################
ioosdif <- SOS(url = "http://opendap.co-ops.nos.noaa.gov/ioos-dif-sos/SOS")
#		method = "GET",
#		verboseOutput = TRUE)
ioosdif_off <- sosOfferings(ioosdif)
sosName(ioosdif_off)

length(ioosdif_off)
# 771

# 
sosResponseFormats(ioosdif)
# TODO test with tab seperated values

#sosObservedProperties(ioosdif_off)
unique(unlist(sosObservedProperties(ioosdif_off)))

# TODO create a map using background data as described in
# http://r-sig-geo.2731867.n2.nabble.com/Reading-a-WMS-layer-as-a-sp-object-td6311972.html


################################################################################
# more data available
#TODO use new data values for mapping example and custom filters
#TODO use KML data from SOS with kmlPlot?

#CO-OPS has expanded its SOS services with an addition of the following 3 new services.
#
#High Low Tide Predictions
#Relative Humidity
#Rain Fall
#
#
#These services are offered for single station and as collections and in the following data formats; CSV, TSV and KML.
#Please note that these services are presently available on the evaluation test site (http://opendap.co-ops.nos.noaa.gov/ioos-dif-sos-test/).
#On January 24th, 2012 at 10:00 am EST, CO-OPS will add these new changes to our operational SOS web site (http://opendap.co-ops.nos.noaa.gov/ioos-dif-sos/).
