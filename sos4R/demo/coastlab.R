# Copyright (C) 2011 by 52 North Initiative for Geospatial Open Source Software GmbH, Contact: info@52north.org
# This program is free software; you can redistribute and/or modify it under the terms of the GNU General Public License version 2 as published by the Free Software Foundation. This program is distributed WITHOUT ANY WARRANTY; even without the implied WARRANTY OF MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program (see gpl-2.0.txt). If not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA or visit the Free Software Foundation web page, http://www.fsf.org.
# Author: Daniel Nuest (daniel.nuest@uni-muenster.de)
# Project: sos4R - visit the project web page, http://www.nordholmen.net/sos4r
library("sos4R")


##############################################################################
# Coastlab SOS - http://www.coastlab.org/
#
# This SOS provides COSYNA data, see the project website at
# http://www.hzg.de/institute/coastal_research/structure/operational_systems/KOK/projects/ICON/index.html
#
# GEOSS description:
# https://geossregistries.info/geosspub/service_details_ns.jsp?serviceId=urn:uuid:b7e0c9d4-9c4a-428c-adb9-fc0efebc9798
#
# Data disclaimer: http://www.coastlab.org/Disclaimer.html
#
# Web Portal: http://kofserver2.hzg.de/codm/
# You can also plot the data here: http://tsdata.hzg.de/index.cgi?seite=plot_form
#
coastlab <- SOS(url = "http://kopc02.gkss.de/sos/sos.py",
		method = SosSupportedConnectionMethods()[["GET"]])

coastlab.off <- sosOfferings(coastlab)
names(coastlab.off)

jade1 <- coastlab.off[["Pile_Jade1"]]
sosObservedProperties(jade1)

jade1.watertemp <- getObservation(sos = coastlab, offering = jade1,
		observedProperty = list("WaterTemperature"), verbose = TRUE,
		eventTime = sosCreateEventTimeList(sosCreateTimePeriod(
						sos = coastlab,
						begin = as.POSIXct(Sys.time() - 3600 * 24), #* 180),
						end = as.POSIXct(Sys.time()))))

# TODO continue implemenation for Coastlab SOS, problem:
# unhandled response document, it contains om:resultDefinition ...
#
#<?xml version="1.0" encoding="UTF-8"?>
#		<om:Observation xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:swe="http://www.opengis.net/swe/0" xmlns:gml="http://www.opengis.net/gml" xmlns:om="http://www.opengis.net/om" xmlns:xlink="http://www.w3.org/1999/xlink" xsi:schemaLocation="http://www.opengis.net/om http://amb25.stccmop.org/schemas/sos/current/sosGetObservation.xsd" gml:id="WaterTemperature">
#		<gml:description>None</gml:description>
#		<gml:name>Pile_Jade1</gml:name>
#		<gml:location>
#		<!-- geometry containing all tuples in this observation -->
#		<gml:Envelope>
#		<gml:lowerCorner srsName="urn:ogc:def:crs:EPSG:6.5:4326">53.516566 8.188217</gml:lowerCorner>
#		<gml:upperCorner srsName="urn:ogc:def:crs:EPSG:6.5:4326">53.516566 8.188217</gml:upperCorner>
#		</gml:Envelope>
#		</gml:location>
#		<!--  Time of response- use TimePeriod for a series of data  -->
#		<!--  or TimeInstant for a single measurement  -->
#		
#		<gml:TimePeriod gml:id="DATA_TIME">
#		<gml:beginPosition>2010-09-19T13:21:38</gml:beginPosition>
#		<gml:endPosition>2011-03-18T12:21:38</gml:endPosition>
#		</gml:TimePeriod>
#		
#		<!-- Procedure -->
#		<om:procedure/>
#		<!-- the property measured -->
#		<om:observedProperty xlink:href="WaterTemperature"/>
#		<!-- Feature Of Interest -->
#		<om:featureOfInterest xlink:href="urn:bodyofwater"/>
#		
#		<!-- Result Structure and Encoding -->
#		<om:resultDefinition>
#		<swe:DataBlockDefinition>
#		<swe:components name="WaterTemperature">
#		<swe:DataRecord>
#		<swe:field name="time">
#		<swe:Time definition="urn:ogc:def:phenomenon:time:iso8601"/>
#		
#		</swe:field>
#		<swe:field name="latitude">
#		<swe:Quantity definition="urn:ogc:def:phenomenon:latitude:wgs84">
#		<swe:uom code="deg"/>
#		</swe:Quantity>
#		</swe:field> 
#		<swe:field name="longitude"> 
#		<swe:Quantity definition="urn:ogc:def:phenomenon:longitude:wgs84">
#		<swe:uom code="deg"/>
#		</swe:Quantity>
#		
#		</swe:field>
#		<swe:field name="depth">
#		<swe:Quantity definition="cf:depth">
#		<swe:uom code="urn:ogc:unit:meter"/>
#		</swe:Quantity>
#		</swe:field>
#		<swe:field name="WaterTemperature">
#		<swe:Quantity definition="WaterTemperature">
#		<swe:uom xlink:href="urn:mm.def:units#degrees C"/>
#		</swe:Quantity>
#		
#		</swe:field>
#		</swe:DataRecord>
#		
#		</swe:components>
#		<swe:encoding>
#		<swe:AsciiBlock tokenSeparator="," blockSeparator="|" decimalSeparator="."/>
#		</swe:encoding>
#		</swe:DataBlockDefinition>
#		</om:resultDefinition>
#		
#		<om:result>2010-09-19T13:30:00Z,53.516566,8.188217,None,14.21|2010-09-19T13:40:00Z,53.516566,8.188217,None,14.21|...</om:result>
#		</om:Observation>

# TODO create parsing function for om:resultDefinition with exchangeability feature!

jade1.watertemp.result <- sosResult(jade1.watertemp)
summary(jade1.watertemp.result)

###################################
# Demo finished, try another one! #
###################################
