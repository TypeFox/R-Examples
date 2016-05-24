################################################################################
# Copyright (C) 2010 by 52 North                                               #
# Initiative for Geospatial Open Source Software GmbH                          #
#                                                                              #
# Contact: Andreas Wytzisk                                                     #
# 52 North Initiative for Geospatial Open Source Software GmbH                 #
# Martin-Luther-King-Weg 24                                                    #
# 48155 Muenster, Germany                                                      #
# info@52north.org                                                             #
#                                                                              #
# This program is free software; you can redistribute and/or modify it under   #
# the terms of the GNU General Public License version 2 as published by the    #
# Free Software Foundation.                                                    #
#                                                                              #
# This program is distributed WITHOUT ANY WARRANTY; even without the implied   #
# WARRANTY OF MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU #
# General Public License for more details.                                     #
#                                                                              #
# You should have received a copy of the GNU General Public License along with #
# this program (see gpl-2.0.txt). If not, write to the Free Software           #
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA or #
# visit the Free Software Foundation web page, http://www.fsf.org.             #
#                                                                              #
# Author: Daniel Nuest (daniel.nuest@uni-muenster.de)                          #
# Created: 2010-06-18                                                          #
# Project: sos4R - visit the project web page, http://www.nordholmen.net/sos4r #
#                                                                              #
################################################################################

################################################################################
# ACTUAL TO STRING AND PRINTING FUNCTIONS

.toString.OwsServiceOperation <- function(x, ...) {
	.s <- paste("Object of class OwsServiceOperation;\n",
			"service: ", x@service, ", version: ", x@version,
			", request: ", x@request)
	return(.s)
}

.print.OwsServiceOperation <- function(x, ...) {
	cat(.toString.OwsServiceOperation(x, ...), "\n")
	invisible(x)
}

.toString.OwsGetCapabilities <- function(x, ...) {
	.s <- paste("Object of class OwsGetCapabilities; service: ",
		x@service, ", request: ", x@request,
		", owsVersion: ", x@owsVersion,
		", acceptVersions: ", toString(x@acceptVersions), sep = " ")

	if (is(x, "OwsGetCapabilities_1.1.0")) {
		.s <- paste(.s, "\n\t sections: ", toString(x@sections),
			", acceptFormats: ", toString(x@acceptFormats),
			", updateSequence: ", x@updateSequence, sep = " ")
			
		if (is(x, "OwsGetCapabilities_2.0.0")) {
			.s <- paste(.s, "\n\tacceptLanguages: ", 
				x@acceptLanguages, sep = " ")
		}
	}
	
	return(.s)
}


.toString.OwsGetCapabilities_1.1.0 <- function(x, ...) {
	return(.toString.OwsGetCapabilities(x, ...), "\n")
}

.print.OwsGetCapabilities <- function(x, ...) {
	cat(.toString.OwsGetCapabilities(x, ...), "\n")
	invisible(x)
}

.print.OwsGetCapabilities_1.1.0 <- function(x, ...) {
	cat(.toString.OwsGetCapabilities(x, ...), "\n")
	invisible(x)
}

.toString.OwsOperationsMetadata <- function(x, ...) {
	.s <- paste("Object of class OwsOperationsMetadata:\n",
			"-- operations:\n", toString(x@operations),
			"-- parameters:\n", toString(x@parameters),
			"-- extendedCapabilities:\n", toString(x@extendedCapabilities))
	return(.s)
}

.print.OwsOperationsMetadata <- function(x, ...) {
	cat(.toString.OwsOperationsMetadata(x, ...), "\n")
	invisible(x)
}

.toString.OwsOperation <- function(x, ...) {
	.s <- paste("Object of class OwsOperation: Name: ", x@name,
			"\n\tParameters (names): ", paste(names(x@parameters)),
			"\n\tDCPs (types): ", paste(names(x@DCPs)),
			"\n\tConstraints: ", toString(x@constraints),
			"\n\tMetadata: ", toString(x@metadata))
	return(.s)
}

.print.OwsOperation <- function(x, ...) {
	cat(.toString.OwsOperation(x, ...), "\n")
	invisible(x)
}

.toString.OwsServiceIdentification <- function(x, ...) {
	.s <- paste("Object of class OwsServiceIdentification:",
			"\n\tServiceType: ", x@serviceType, "; serviceTypeVersion(s): ",
			paste(x@serviceTypeVersion, collapse = ", "),
			"\n\ttitle(s): ", paste(x@title, collapse = "; "),
			# optional:
			"\n\tProfile(s): ", toString(x@profile),
			"\n\tAbstract(s): ", toString(x@abstract),
			"\n\tKeywords(s): ", toString(x@keywords),
			"\n\tAccessConstraints(s): ", paste(x@accessConstraints,
					collapse = "; "))
	return(.s)
}

.print.OwsServiceIdentification <- function(x, ...) {
	cat(.toString.OwsServiceIdentification(x, ...), "\n")
	invisible(x)
}

.toString.OwsServiceProvider <- function(x, ...) {
	.s <- paste("Object of class OwsServiceProvider:\n\tProvider name: ",
			x@providerName,	"; providerSite: ", x@providerSite,
			"\n\tService contact:  (unparsed XML, see @serviceContact for details)")
	return(.s)
}

.print.OwsServiceProvider <- function(x, ...) {
	cat(.toString.OwsServiceProvider(x, ...), "\n")
	invisible(x)
}

.toString.OwsContents <- function(x, ...) {
	return("Object of class OwsContents (wraps unparsed XML, see @xml for details).\n")
}

.print.OwsContents <- function(x, ...) {
	cat(.toString.OwsContents(x, ...), "\n")
	invisible(x)
}

.toString.OwsCapabilities <- function(x, ...) {
	.s <- paste("Object of class ", class(x), " -- version: ", x@version,
			", owsVersion: ", x@owsVersion, ", updateSequence: ", 
			x@updateSequence)
	return(.s)
}

.print.OwsCapabilities <- function(x, ...) {
	cat(.toString.OwsCapabilities(x, ...), "\n")
	invisible(x)
}

.toString.OwsExceptionReport <- function(x, ...) {
	.s <- paste(
			"Object of class OwsExceptionReport",
			"; version: ",
			x@version,
			"; lang: ",
			x@lang,
			paste(";\n", length(x@exceptions),
					"exception(s) (code @ locator : text):"),
			sep = "")
	for (e in x@exceptions) {
		.e <- paste("  ", e@exceptionCode, " @ ", e@locator, " :\n\t",
				e@exceptionText, sep = "")
		.s <- paste(.s, .e, sep = "\n")
	}
	return(paste(.s, "\n"))
}

.print.OwsExceptionReport <- function(x, ...) {
	cat(.toString.OwsExceptionReport(x, ...), "\n")
	invisible(x)
}

.toString.OwsException <- function(x, ...) {
	.s <- paste("Object of class OwsException; exception code: ",
			x@exceptionCode,
			", locator: ",
			x@locator,
			"\nException text(s):\n\t",
			paste(x@exceptionText, sep = "\n\t"))
	return(.s)
}

.print.OwsException <- function(x, ...) {
	cat(.toString.OwsException, "\n")
	invisible(x)
}

.toString.OwsRange <- function(x, ...) {
	.s <- paste("Object of class OwsRange; spacing: ",
			x@spacing, ", rangeClosure: ",
			x@rangeClosure,
			"\nFROM ", x@minimumValue, " TO ",x@maximumValue, sep = "")
	return(.s)
}

.print.OwsRange <- function(x, ...) {
	cat(.toString.OwsRange(x, ...), "\n")
	invisible(x)
}

.toString.SOS <- function(x, ...) {
	.s <- paste("Object of class SOS -- version: ", x@version,
			"\n\tCapabilities: ", toString(x@capabilities))
	return(.s)
}

.print.SOS <- function(x, ...) {
	cat(.toString.SOS(x, ...), "\n")
	invisible(x)
}

.toString.SOS_1.0.0 <- function(x, ...) {
	.s <- paste("Object of class SOS_1.0.0 [",
			x@method,
			", ",
			x@url,
			", ",
			sosTitle(x),
			"]", sep = "")
	return(.s)
}

.print.SOS_1.0.0 <- function(x, ...) {
	cat(.toString.SOS_1.0.0(x, ...), "\n")
	invisible(x)
}

.toString.SosFilter_Capabilities <- function(x, ...) {
	.s <- paste("Object of class SosFilter_Capabilities;",
			"\n\tSpatial_Capabilities:\t",
			toString(x@spatial[[ogcGeometryOperandsName]]),
			";",
			toString(x@spatial[[ogcSpatialOperatorsName]]),
			"\n\tTemporal_Capablities:\t",
			toString(x@temporal[[ogcTemporalOperandsName]]),
			";",
			toString(x@spatial[[ogcTemporalOperatorsName]]),
			"\n\tScalar_Capablities:\t\t",
			toString(x@scalar),
			"\n\tId_Capabilities\t\t\t",
			toString(x@id))
	return(.s)
}

.print.SosFilter_Capabilities <- function(x, ...) {
	cat(.toString.SosFilter_Capabilities(x, ...), "\n")
	invisible(x)
}

.toString.SosObservationOffering <- function(x, ...) {
	.s <- paste("Object of class SosObservationOffering; ",
			"id: ", x@id, ", name: ", x@name,
			"\n\ttime: ", .addTabIndent(toString(x@time)),
			"\n\tprocedure(s): ", toString(paste(x@procedure)),
			"\n\tobservedProperty(s): ", toString(paste(x@observedProperty)),
			"\n\tfeature(s)OfInterest: ", toString(paste(x@featureOfInterest)),
			"\n\tresponseFormat(s): ",  toString(paste(x@responseFormat)),
			", responseMode(s): ",  toString(paste(x@responseMode)),
			"\n\tintendedApplication: ", toString(x@intendedApplication),
			"\n\tresultModel(s): ",  toString(paste(x@resultModel)),
			 "\n\tboundedBy: ",  toString(paste(x@boundedBy)))
	return(.s)
}

.print.SosObservationOffering <- function(x, ...) {
	cat(.toString.SosObservationOffering(x, ...), "\n")
	invisible(x)
}

.toString.SosContents <- function(x, ...) {
	.s <- paste("Object of class SosContents with observation offerings (names):\n\t",
			toString(paste(names(x@observationOfferings))))
	return(.s)
}

.print.SosContents <- function(x, ...) {
	cat(.toString.SosContents(x, ...), "\n")
	invisible(x)
}

.toString.SosEventTime <- function(x, ...) {
	.s <- paste("Object of class SosEventTime:\n\t",
			class(x@temporalOps),": ",
			toString(x@temporalOps@time), sep = "")
	return(.s)
}

.print.SosEventTime <- function(x, ...) {
	cat(.toString.SosEventTime(x, ...), "\n")
	invisible(x)
}

.toString.SosEventTimeLatest <- function(x, ...) {
	.s <- paste("Object of class SosEventTimeLatest; temporalOps value:",
			x@temporalOps)
	return(.s)
}

.print.SosEventTimeLatest <- function(x, ...) {
	cat(.toString.SosEventTimeLatest(x, ...), "\n")
	invisible(x)
}

.toString.SosFeatureOfInterest <- function(x, ...) {
	.s <- paste("Object of class SosFeatureOfInterest",
			";\n\tobjectIDs: ",
			toString(x@objectIDs),
			";\n\tspatialOps: ",
			toString(x@spatialOps),
			sep = "")
	return(.s)
}

.print.SosFeatureOfInterest <- function(x, ...) {
	cat(.toString.SosFeatureOfInterest(x, ...), "\n")
	invisible(x)
}

.toString.SensorML <- function(x, ...) {
	.s <- ("Object of class SensorML (see @xml for full document).\n")
	.s <- paste(.s, "\tID: ", sosId(x),
			"\n\tname:", sosName(x),
			"\n\tdescription:", sosAbstract(x),
			"\n\tcoords:", toString(sosCoordinates(x)),
			"\n\tboundedBy:", toString(sosBoundedBy(x)))
	return(.s)
}

.print.SensorML <- function(x, ...) {
	cat(.toString.SensorML(x, ...), "\n")
	invisible(x)
}

.toString.SosGetObservation <- function(x, ...) {
	.s <- paste("Object of class SosGetObservation: ",
			"service: ",
			x@service,
			", version: ",
			x@version,
			", offering: ",
			x@offering,
			"\nobservered property: ",
			x@observedProperty,
			"\nresponseFormat(s): ",
			x@responseFormat,
			", responseMode(s): ",
			toString(paste(x@responseMode)),
			# optionals:
			"\nprocedure(s)",
			toString(paste(x@procedure)),
			"\n\tfeature(s) of interest",
			toString(x@featureOfInterest),
			"\n\tevent time: ",
			toString(x@eventTime),
			"\n\tresult: ",
			class(x@result),
			"\nsrsName: ",
			x@srsName,
			"\nresultModel(s): ",
			x@resultModel)
	return(.s)
}

.print.SosGetObservation <- function(x, ...) {
	cat(.toString.SosGetObservation(x, ...), "\n")
	invisible(x)
}

.toString.SosGetObservationById <- function(x, ...) {
	.s <- paste("Object of class SosGetObservationById: ",
			"service: ",
			x@service,
			", version: ",
			x@version,
			"\nObsvervation ID: ",
			x@observationId,
			"\nResponseFormat(s): ",
			x@responseFormat,
			", responseMode(s): ",
			toString(paste(x@responseMode)),
			# optionals:
			", srsName: ",
			x@srsName,
			", resultModel(s): ",
			x@resultModel)
	return(.s)
}

.print.SosGetObservationById <- function(x, ...) {
	cat(.toString.SosGetObservationById(x, ...), "\n")
	invisible(x)
}


.toString.SosDescribeSensor <- function(x, ...) {
	.s <- paste("Object of class SosDescribeSensor: ",
			"service: ",
			x@service,
			", version: ",
			x@version,
			", outputFormat: ",
			x@outputFormat,
			"\nProcedure: ",
			x@procedure)
	return(.s)
}

.print.SosDescribeSensor <- function(x, ...) {
	cat(.toString.SosDescribeSensor(x, ...), "\n")
	invisible(x)
}

.toString.OmMeasurement <- function(x, ...) {
	.s <- paste(
			"Object of class OmMeasurement",
			", procedure ",
			toString(x@procedure),
			", observedProperty: ",
			toString(x@observedProperty),
			";\n\tfeatureOfInterest: ",
			toString(x@featureOfInterest),
			";\n\tsamplingTime: ",
			toString(x@samplingTime),
			";\n\tresult: ",
			toString(x@result),
			sep = "")
	return(.s)
}

.print.OmMeasurement <- function(x, ...) {
	cat(.toString.OmMeasurement(x, ...), "\n")
	invisible(x)
}

.toString.OmObservationCollection <- function(x, ...) {
	.s <- paste("Object of class OmObservationCollection with",
			length(x), "members.")
	return(.s)
}

.print.OmObservationCollection <- function(x, ...) {
	cat(.toString.OmObservationCollection(x, ...), "\n")
	invisible(x)
}

.toString.OmObservation <- function(x, ...) {
	if(is.null(sosObservedProperties(x)))
		.obsProp <- NULL
	else .obsProp <- toString(sosObservedProperties(x))
	
	.s <- paste("Object of class OmObservation;\n\tprocedure: ",
			toString(x@procedure),
			"\n\tobservedProperty: ", .obsProp,
			"\n\tfoi: ", toString(sosFeatureIds(x)),
			"\n\tsamplingTime: ", .addTabIndent(toString(x@samplingTime)),
			"\n\tresult dimensions: ", toString(dim(x@result)),
			sep = "")
	return(.s)
}

.print.OmObservation <- function(x, ...) {
	cat(.toString.OmObservation(x, ...), "\n")
	invisible(x)
}

.toString.OmObservationProperty <- function(x, ...) {
	.s <- paste("Object of class OmObservationProperty",
			"; href: ",
			x@href,
			"; observation: ",
			toString(x@obs),
			sep = "")
	return(.s)
}

.print.OmObservationProperty <- function(x, ...) {
	cat(.toString.OmObservationProperty(x, ...), "\n")
	invisible(x)
}

.toString.GmlMeasure <- function(x, ...) {
	.s <- paste(
			"Object of class GmlMeasure",
			"; value: ",
			x@value,
			"; uom: ",
			x@uom,
			sep = "")
	return(.s)
}

.print.GmlMeasure <- function(x, ...) {
	cat(.toString.GmlMeasure(x, ...), "\n")
	invisible(x)
}

.toString.SwePhenomenonProperty <- function(x, ...) {
	.s <- paste("Object of class SwePhenomenonProperty",
			"; href: ",
			x@href,
			"; phenomenon: ",
			toString(x@phenomenon),
			sep = "")
	return(.s)
}

.print.SwePhenomenonProperty <- function(x, ...) {
	cat(.toString.SwePhenomenonProperty(x, ...), "\n")
	invisible(x)
}

.toString.SwePhenomenon <- function(x, ...) {
	.s <- paste("Object of class SwePhenomenon",
		"; id: ",
		x@id,
		"; name: ",
		x@name,
		"; description: ",
		x@description,
		sep = "")
	invisible(x)
}

.print.SwePhenomenon <- function(x, ...) {
	cat(.toString.SwePhenomenon(x, ...), "\n")
	invisible(x)
}

.toString.SweCompositePhenomenon <- function(x, ...) {
	.s <- paste("Object of class SweCompositePhenomenon",
		"; id: ",
		x@id,
		"; name: ",
		x@name,
		"; description: ",
		x@description,
		"; dimension: ",
		x@dimension,
		"; base: ",
		toString(x@base),
		";\ncomponents:\t",
		sapply(sapply(x@components, toString), paste, "\t\t\t"),
		sep = "")
	return(.s)
}

.print.SweCompositePhenomenon <- function(x, ...) {
	cat(.toString.SweCompositePhenomenon(x, ...), "\n")
	invisible(x)
}

.toString.SweTextBlock <- function(x, ...) {
	.s <- paste("Object of class SweTextBlock",
			" '",
			x@tokenSeparator,
			" ",
			x@blockSeparator,
			" ",
			x@decimalSeparator,
			"'; id: ",
			x@id,
			sep = "")
	return(.s)
}

.print.SweTextBlock <- function(x, ...) {
	cat(.toString.SweTextBlock(x, ...), "\n")
	invisible(x)
}

.toString.GmlTimePosition <- function(x, ...) {
	.s <- paste("GmlTimePosition [",
			" time: ", x@time, sep = "")
	
	if(!is.na(x@frame)) {
		.s <- paste(.s, "; frame: ", x@frame, sep = "")
	}
	if(!is.na(x@frame)) {
		.s <- paste(.s, "; frame: ", x@frame, sep = "")
	}
	if(!is.na(x@calendarEraName)) {
		.s <- paste(.s, "; calendarEraName: ", x@calendarEraName, sep = "")
	}
	if(!is.na(x@indeterminatePosition)) {
		.s <- paste(.s, "; indeterminatePosition: ", x@indeterminatePosition,
				sep = "")
	}
	
	.s <- paste(.s, "]")
	return(.s)
}

.print.GmlTimePosition <- function(x, ...) {
	cat(.toString.GmlTimePosition(x, ...), "\n")
	invisible(x)
}

.toString.GmlTimeInstant <- function(x, ...) {
	.s <- paste(
			#"Object of class GmlTimeInstant",
			#"; timePosition:",
			toString(x@timePosition))
	return(.s)
}

.print.GmlTimeInstant <- function(x, ...) {
	cat(.toString.GmlTimeInstant(x, ...))
	invisible(x)
}

.toString.GmlTimeInstantProperty <- function(x, ...) {
	.s <- paste("Object of class GmlTimeInstantProperty",
			"; href: ",
			x@href,
			"; time: ",
			toString(x@time),
			sep = "")
	return(.s)
}

.print.GmlTimeInstantProperty <- function(x, ...) {
	cat(.toString.GmlTimeInstantProperty(x, ...), "\n")
	invisible(x)
}

.toString.GmlTimeInterval <- function(x, ...) {
	.s <- paste("Object of class GmlTimeInterval",
			"; interval: ",
			x@interval,
			"; unit: ",
			x@unit,
			"; radix: ",
			x@radix,
			"; factor: ",
			x@factor,
			sep = "")
	return(.s)
}

.print.GmlTimeInterval <- function(x, ...) {
	cat(.toString.GmlTimeInterval(x, ...), "\n")
	invisible(x)
}

.toString.GmlTimePeriod <- function(x, ...) {
	.s <- ""

	if(!is.na(x@duration)) {
		.s <- paste(.s, "; duration: ", x@duration)
	}
	if(!is.null(x@timeInterval)) {
		.s <- paste(.s, ", timeInterval: ", toString(x@timeInterval), ";",
				sep = "")
	}
	
	if(!is.null(x@begin) && !is.null(x@end)) {
		.s <- paste(toString(x@begin), "\n\t--> ", toString(x@end), sep = "")
	}
	else {
		.s <- paste(toString(x@beginPosition), "\n\t--> ",
				toString(x@endPosition), sep = "")
	}
	
	.s <- paste("GmlTimePeriod: [", .s, "]")
	return(.s)
}

.print.GmlTimePeriod <- function(x, ...) {
	cat(.toString.GmlTimePeriod(x, ...), "\n")
	invisible(x)
}

.toString.GmlFeatureProperty <- function(x, ...) {
	.s <- paste("Object of class GmlFeatureProperty",
			", href: ",
			x@href,
			", feature: ",
			toString(x@feature),
			sep = "")
	return(.s)
}

.print.GmlFeatureProperty <- function(x, ...) {
	cat(.toString.GmlFeatureProperty(x, ...), "\n")
	invisible(x)
}

.toString.GmlFeatureCollection <- function(x, ...) {
	.s <- paste("Object of class GmlFeatureCollection",
			"; id: ", x@id, ";\n\t",
			length(x@featureMembers), " featureMember(s): ",
			toString(x@featureMembers),
			sep = "")
	return(.s)
}

.print.GmlFeatureCollection <- function(x, ...) {
	cat(.toString.GmlFeatureCollection(x, ...), "\n")
	invisible(x)
}

.toString.GmlDirectPosition <- function(x, ...) {
	.s <- paste("Object of class GmlDirectPosition",
			"; pos: ",
			x@pos,
			"; srsName: ",
			x@srsName,
			", srsDimension: ",
			x@srsDimension,
			", srsAxisLabels: ",
			x@axisLabels,
			", uomLabels: ",
			x@uomLabels,
			sep = "")
	return(.s)
}

.print.GmlDirectPosition <- function(x, ...) {
	cat(.toString.GmlDirectPosition(x, ...), "\n")
	invisible(x)
}

.toString.GmlPoint <- function(x, ...) {
	.s <- paste("Object of class GmlPoint",
			"; pos: ",
			toString(x@pos),
			";\nsrsName: ",
			x@srsName,
			", srsDimension: ",
			x@srsDimension,
			", srsDimension: ",
			x@axisLabels,
			", uomLabels: ",
			x@uomLabels,
			sep = "")
	return(.s)
}

.print.GmlPoint <- function(x, ...) {
	cat(.toString.GmlPoint(x, ...), "\n")
	invisible(x)
}

.toString.GmlPointProperty <- function(x, ...) {
	.s <- paste("Object of class GmlPointProperty",
			"; href: ",
			x@href,
			"; point: ",
			toString(x@point),
			sep = "")
	return(.s)
}

.print.GmlPointProperty <- function(x, ...) {
	cat(.toString.GmlPointProperty(x, ...), "\n")
	invisible(x)
}

.tempOpToString <- function(obj) {
	.s <- paste("propertyName:", obj@propertyName,
			"time:", toString(obj@time))
	return(.s)
}

.toString.GmlGeometry <- function(x, ...) {
	.s <- paste("Object of class GmlGeometry",
			"; id: ",
			x@id,
			";\nsrsName: ",
			x@srsName,
			", srsDimension: ",
			x@srsDimension,
			", srsDimension: ",
			x@axisLabels,
			", uomLabels: ",
			x@uomLabels,
			sep = "")
	return(.s)
}

.print.GmlGeometry <- function(x, ...) {
	cat(.toString.GmlGeometry(x, ...), "\n")
	invisible(x)
}

.toString.GmlEnvelope <- function(x, ...) {
	.s <- paste("Object of class GmlEnvelope",
			"; srsName: ",
			x@srsName,
			", srsDimension: ",
			x@srsDimension,
			", srsDimension: ",
			x@axisLabels,
			", uomLabels: ",
			x@uomLabels,
			";\n\tlowerCorner: ",
			toString(x@lowerCorner),
			";\n\tupperCorner: ",
			toString(x@upperCorner),
			sep = "")
	return(.s)
}

.print.GmlEnvelope <- function(x, ...) {
	cat(.toString.GmlEnvelope(x, ...), "\n")
	invisible(x)
}

.toString.TM_After <- function(x, ...) {
	.s <- paste("Object of class TM_After;",
			.tempOpToString(x))
	return(.s)
}

.print.TM_After <- function(x, ...) {
	cat(.toString.TM_After(x, ...), "\n")
	invisible(x)
}

.toString.TM_Before <- function(x, ...) {
	.s <- paste("Object of class TM_Before;",
			.tempOpToString(x))
	return(.s)
}

.print.TM_Before <- function(x, ...) {
	cat(.toString.TM_Before(x, ...), "\n")
	invisible(x)
}

.toString.TM_During <- function(x, ...) {
	.s <- paste("Object of class TM_During;",
			.tempOpToString(x))
	return(.s)
}

.print.TM_During <- function(x, ...) {
	cat(.toString.TM_During(x, ...), "\n")
	invisible(x)
}

.toString.TM_Equals <- function(x, ...) {
	.s <- paste("Object of class TM_Equals;",
			.tempOpToString(x))
	return(.s)
}


.print.TM_Equals <- function(x, ...) {
	cat(.toString.TM_Equals(x, ...), "\n")
	invisible(x)
}

.toString.OgcBBOX <- function(x, ...) {
	.s <- paste("Object of class OgcBBOX; propertyName: ",
			x@propertyName,
			"; envelope: ",
			toString(x@envelope),
			sep = "")
	return(.s)
}

.print.OgcBBOX <- function(x, ...) {
	cat(.toString.OgcBBOX(x, ...), "\n")
	invisible(x)
}

.binSpatOpToString <- function(x, ...) {
	.s <- paste("propertyName:",
			x@propertyName,
			";\n\tgeometry: ",
			toString(x@geometry),
			";\n\tenvelope: ",
			toString(x@envelope),
			sep = "")
	return(.s)
}

.toString.OgcContains <- function(x, ...) {
	.s <- paste("Object of class OgcContains;",
			.binSpatOpToString(x))
	return(.s)
}

.print.OgcContains <- function(x, ...) {
	cat(.toString.OgcContains(x, ...), "\n")
	invisible(x)
}

.toString.OgcIntersects <- function(x, ...) {
	.s <- paste("Object of class OgcIntersects;",
			.binSpatOpToString(x))
	return(.s)
}

.print.OgcIntersects <- function(x, ...) {
	cat(.toString.OgcIntersects(x, ...), "\n")
	invisible(x)
}

.toString.OgcOverlaps <- function(x, ...) {
	.s <- paste("Object of class OgcOverlaps;",
			.binSpatOpToString(x))
	return(.s)
}

.print.OgcOverlaps <- function(x, ...) {
	cat(.toString.OgcOverlaps(x, ...), "\n")
	invisible(x)
}

.toString.SaSamplingPoint <- function(x, ...) {
	.s <- paste("Object of class SaSamplingPoint",
			"; id: ",
			x@id,
			"; position: ",
			toString(x@position),
			", relatedObservation: ",
			toString(x@relatedObservation),
			", relatedSamplingFeature: ",
			toString(x@relatedSamplingFeature),
			", surveyDetails: ",
			toString(x@surveyDetails),
			";\n\tsampledFeatures: ",
			toString(x@sampledFeatures),
			sep = "")
	return(.s)
}

.print.SaSamplingPoint <- function(x, ...) {
	cat(.toString.SaSamplingPoint(x, ...), "\n")
	invisible(x)
}

.toString.SaSamplingSurface <- function(x, ...) {
	.s <- paste("Object of class SaSamplingSurface",
			"; id: ",
			x@id,
			"; shape: ",
			toString(x@shape),
			", relatedObservation: ",
			toString(x@relatedObservation),
			", relatedSamplingFeature: ",
			toString(x@relatedSamplingFeature),
			", surveyDetails: ",
			toString(x@surveyDetails),
			", position: ",
			toString(x@position),
			";\n\tsampledFeatures: ",
			toString(x@sampledFeatures),
			sep = "")
	return(.s)
}

.print.SaSamplingSurface <- function(x, ...) {
	cat(.toString.SaSamplingSurface(x, ...), "\n")
	invisible(x)
}

################################################################################
# PRINT FUNCTIONS
setMethod("print", "OwsServiceOperation", function(x, ...) .print.OwsServiceOperation(x, ...))
setMethod("print", "OwsGetCapabilities", function(x, ...) .print.OwsGetCapabilities(x, ...))
setMethod("print", "OwsGetCapabilities_1.1.0", function(x, ...) .print.OwsGetCapabilities(x, ...))
setMethod("print", "OwsGetCapabilities_2.0.0", function(x, ...) .print.OwsGetCapabilities(x, ...))
setMethod("print", "OwsOperationsMetadata", function(x, ...) .print.OwsOperationsMetadata(x, ...))
setMethod("print", "OwsOperation", function(x, ...) .print.OwsOperation(x, ...))
setMethod("print", "OwsServiceIdentification", function(x, ...) .print.OwsServiceIdentification(x, ...))
setMethod("print", "OwsServiceProvider", function(x, ...) .print.OwsServiceProvider(x, ...))
setMethod("print", "OwsContents", function(x, ...) .print.OwsContents(x, ...))
setMethod("print", "OwsCapabilities", function(x, ...) .print.OwsCapabilities(x, ...))
setMethod("print", "OwsCapabilities_1.1.0", function(x, ...) .print.OwsCapabilities(x, ...))
setMethod("print", "OwsCapabilities_2.0.0", function(x, ...) .print.OwsCapabilities(x, ...))
setMethod("print", "OwsExceptionReport", function(x, ...) .print.OwsExceptionReport(x, ...))
setMethod("print", "OwsException", function(x, ...) .print.OwsException(x, ...))
setMethod("print", "OwsRange", function(x, ...) .print.OwsRange(x, ...))
setMethod("print", "SOS", function(x, ...) .print.SOS(x, ...))
setMethod("print", "SOS_1.0.0", function(x, ...) .print.SOS_1.0.0(x, ...))
setMethod("print", "SosFilter_Capabilities", function(x, ...) .print.SosFilter_Capabilities(x, ...))
setMethod("print", "SosObservationOffering", function(x, ...) .print.SosObservationOffering(x, ...))
setMethod("print", "SosContents", function(x, ...) .print.SosContents(x, ...))
setMethod("print", "SosEventTime", function(x, ...) .print.SosEventTime(x, ...))
setMethod("print", "SosFeatureOfInterest", function(x, ...) .print.SosFeatureOfInterest(x, ...))
setMethod("print", "SensorML", function(x, ...) .print.SensorML(x, ...))
setMethod("print", "SosGetObservation", function(x, ...) .print.SosGetObservation(x, ...))
setMethod("print", "SosGetObservationById", function(x, ...) .print.SosGetObservationById(x, ...))
setMethod("print", "SosDescribeSensor", function(x, ...) .print.SosDescribeSensor(x, ...))
setMethod("print", "SaSamplingPoint", function(x, ...) .print.SaSamplingPoint(x, ...))
setMethod("print", "SaSamplingSurface", function(x, ...) .print.SaSamplingSurface(x, ...))
setMethod("print", "SwePhenomenon", function(x, ...) .print.SwePhenomenon(x, ...))
setMethod("print", "SwePhenomenonProperty", function(x, ...) .print.SwePhenomenonProperty(x, ...))
setMethod("print", "SweCompositePhenomenon", function(x, ...) .print.SweCompositePhenomenon(x, ...))
setMethod("print", "SweTextBlock", function(x, ...) .print.SweTextBlock(x, ...))
setMethod("print", "OmObservationCollection", function(x, ...) .print.OmObservationCollection(x, ...))
setMethod("print", "OmObservation", function(x, ...) .print.OmObservation(x, ...))
setMethod("print", "OmObservationProperty", function(x, ...) .print.OmObservationProperty(x, ...))
setMethod("print", "GmlMeasure", function(x, ...) .print.GmlMeasure(x, ...))
setMethod("print", "OmMeasurement", function(x, ...) .print.OmMeasurement(x, ...))
setMethod("print", "GmlTimePosition", function(x, ...) .print.GmlTimePosition(x, ...))
setMethod("print", "GmlTimeInstant", function(x, ...) .print.GmlTimeInstant(x, ...))
setMethod("print", "GmlTimeInstantProperty", function(x, ...) .print.GmlTimeInstantProperty(x, ...))
setMethod("print", "GmlTimeInterval", function(x, ...) .print.GmlTimeInterval(x, ...))
setMethod("print", "GmlTimePeriod", function(x, ...) .print.GmlTimePeriod(x, ...))
setMethod("print", "GmlFeatureProperty", function(x, ...) .print.GmlFeatureProperty(x, ...))
setMethod("print", "GmlFeatureCollection", function(x, ...) .print.GmlFeatureCollection(x, ...))
setMethod("print", "GmlDirectPosition", function(x, ...) .print.GmlDirectPosition(x, ...))
setMethod("print", "GmlPoint", function(x, ...) .print.GmlPoint(x, ...))
setMethod("print", "GmlPointProperty", function(x, ...) .print.GmlPointProperty(x, ...))
setMethod("print", "GmlGeometry", function(x, ...) .print.GmlGeometry(x, ...))
setMethod("print", "GmlEnvelope", function(x, ...) .print.GmlEnvelope(x, ...))
setMethod("print", "TM_After", function(x, ...) .print.TM_After(x, ...))
setMethod("print", "TM_Before", function(x, ...) .print.TM_Before(x, ...))
setMethod("print", "TM_During", function(x, ...) .print.TM_During(x, ...))
setMethod("print", "TM_Equals", function(x, ...) .print.TM_Equals(x, ...))
setMethod("print", "OgcBBOX", function(x, ...) .print.OgcBBOX(x, ...))
setMethod("print", "OgcContains", function(x, ...) .print.OgcContains(x, ...))
setMethod("print", "OgcIntersects", function(x, ...) .print.OgcIntersects(x, ...))
setMethod("print", "OgcOverlaps", function(x, ...) .print.OgcOverlaps(x, ...))

################################################################################
# TO STRING FUNCTIONS
setMethod("toString", "OwsServiceOperation", function(x, ...) .toString.OwsServiceOperation(x, ...))
setMethod("toString", "OwsGetCapabilities", function(x, ...) .toString.OwsGetCapabilities(x, ...))
setMethod("toString", "OwsGetCapabilities_1.1.0", function(x, ...) .toString.OwsGetCapabilities(x, ...))
setMethod("toString", "OwsGetCapabilities_2.0.0", function(x, ...) .toString.OwsGetCapabilities(x, ...))
setMethod("toString", "OwsOperationsMetadata", function(x, ...) .toString.OwsOperationsMetadata(x, ...))
setMethod("toString", "OwsOperation", function(x, ...) .toString.OwsOperation(x, ...))
setMethod("toString", "OwsServiceIdentification", function(x, ...) .toString.OwsServiceIdentification(x, ...))
setMethod("toString", "OwsServiceProvider", function(x, ...) .toString.OwsServiceProvider(x, ...))
setMethod("toString", "OwsContents", function(x, ...) .toString.OwsContents(x, ...))
setMethod("toString", "OwsCapabilities", function(x, ...) .toString.OwsCapabilities(x, ...))
setMethod("toString", "OwsCapabilities_1.1.0", function(x, ...) .toString.OwsCapabilities(x, ...))
setMethod("toString", "OwsCapabilities_2.0.0", function(x, ...) .toString.OwsCapabilities(x, ...))
setMethod("toString", "OwsExceptionReport", function(x, ...) .toString.OwsExceptionReport(x, ...))
setMethod("toString", "OwsException", function(x, ...) .toString.OwsException(x, ...))
setMethod("toString", "OwsRange", function(x, ...) .toString.OwsRange(x, ...))
setMethod("toString", "SOS", function(x, ...) .toString.SOS(x, ...))
setMethod("toString", "SOS_1.0.0", function(x, ...) .toString.SOS_1.0.0(x, ...))
setMethod("toString", "SosFilter_Capabilities", function(x, ...) .toString.SosFilter_Capabilities(x, ...))
setMethod("toString", "SosObservationOffering", function(x, ...) .toString.SosObservationOffering(x, ...))
setMethod("toString", "SosContents", function(x, ...) .toString.SosContents(x, ...))
setMethod("toString", "SosEventTime", function(x, ...) .toString.SosEventTime(x, ...))
setMethod("toString", "SosFeatureOfInterest", function(x, ...) .toString.SosFeatureOfInterest(x, ...))
setMethod("toString", "SensorML", function(x, ...) .toString.SensorML(x, ...))
setMethod("toString", "SosGetObservation", function(x, ...) .toString.SosGetObservation(x, ...))
setMethod("toString", "SosGetObservationById", function(x, ...) .toString.SosGetObservationById(x, ...))
setMethod("toString", "SosDescribeSensor", function(x, ...) .toString.SosDescribeSensor(x, ...))
setMethod("toString", "SaSamplingPoint", function(x, ...) .toString.SaSamplingPoint(x, ...))
setMethod("toString", "SaSamplingSurface", function(x, ...) .toString.SaSamplingSurface(x, ...))
setMethod("toString", "SwePhenomenon", function(x, ...) .toString.SwePhenomenon(x, ...))
setMethod("toString", "SwePhenomenonProperty", function(x, ...) .toString.SwePhenomenonProperty(x, ...))
setMethod("toString", "SweCompositePhenomenon", function(x, ...) .toString.SweCompositePhenomenon(x, ...))
setMethod("toString", "SweTextBlock", function(x, ...) .toString.SweTextBlock(x, ...))
setMethod("toString", "OmObservationCollection", function(x, ...) .toString.OmObservationCollection(x, ...))
setMethod("toString", "OmObservation", function(x, ...) .toString.OmObservation(x, ...))
setMethod("toString", "OmObservationProperty", function(x, ...) .toString.OmObservationProperty(x, ...))
setMethod("toString", "GmlMeasure", function(x, ...) .toString.GmlMeasure(x, ...))
setMethod("toString", "OmMeasurement", function(x, ...) .toString.OmMeasurement(x, ...))
setMethod("toString", "GmlTimePosition", function(x, ...) .toString.GmlTimePosition(x, ...))
setMethod("toString", "GmlTimeInstant", function(x, ...) .toString.GmlTimeInstant(x, ...))
setMethod("toString", "GmlTimeInstantProperty", function(x, ...) .toString.GmlTimeInstantProperty(x, ...))
setMethod("toString", "GmlTimeInterval", function(x, ...) .toString.GmlTimeInterval(x, ...))
setMethod("toString", "GmlTimePeriod", function(x, ...) .toString.GmlTimePeriod(x, ...))
setMethod("toString", "GmlFeatureProperty", function(x, ...) .toString.GmlFeatureProperty(x, ...))
setMethod("toString", "GmlFeatureCollection", function(x, ...) .toString.GmlFeatureCollection(x, ...))
setMethod("toString", "GmlDirectPosition", function(x, ...) .toString.GmlDirectPosition(x, ...))
setMethod("toString", "GmlPoint", function(x, ...) .toString.GmlPoint(x, ...))
setMethod("toString", "GmlPointProperty", function(x, ...) .toString.GmlPointProperty(x, ...))
setMethod("toString", "GmlGeometry", function(x, ...) .toString.GmlGeometry(x, ...))
setMethod("toString", "GmlEnvelope", function(x, ...) .toString.GmlEnvelope(x, ...))
setMethod("toString", "TM_After", function(x, ...) .toString.TM_After(x, ...))
setMethod("toString", "TM_Before", function(x, ...) .toString.TM_Before(x, ...))
setMethod("toString", "TM_During", function(x, ...) .toString.TM_During(x, ...))
setMethod("toString", "TM_Equals", function(x, ...) .toString.TM_Equals(x, ...))
setMethod("toString", "OgcBBOX", function(x, ...) .toString.OgcBBOX(x, ...))
setMethod("toString", "OgcContains", function(x, ...) .toString.OgcContains(x, ...))
setMethod("toString", "OgcIntersects", function(x, ...) .toString.OgcIntersects(x, ...))
setMethod("toString", "OgcOverlaps", function(x, ...) .toString.OgcOverlaps(x, ...))

################################################################################
# SHOW FUNCTIONS
setMethod("show", "OwsServiceOperation", function(object) .print.OwsServiceOperation(object))
setMethod("show", "OwsGetCapabilities", function(object) .print.OwsGetCapabilities(object))
setMethod("show", "OwsGetCapabilities_1.1.0", function(object) .print.OwsGetCapabilities(object))
setMethod("show", "OwsGetCapabilities_2.0.0", function(object) .print.OwsGetCapabilities(object))
setMethod("show", "OwsOperationsMetadata", function(object) .print.OwsOperationsMetadata(object))
setMethod("show", "OwsOperation", function(object) .print.OwsOperation(object))
setMethod("show", "OwsServiceIdentification", function(object) .print.OwsServiceIdentification(object))
setMethod("show", "OwsServiceProvider", function(object) .print.OwsServiceProvider(object))
setMethod("show", "OwsContents", function(object) .print.OwsContents(object))
setMethod("show", "OwsCapabilities", function(object) .print.OwsCapabilities(object))
setMethod("show", "OwsCapabilities_1.1.0", function(object) .print.OwsCapabilities(object))
setMethod("show", "OwsCapabilities_2.0.0", function(object) .print.OwsCapabilities(object))
setMethod("show", "OwsExceptionReport", function(object) .print.OwsExceptionReport(object))
setMethod("show", "OwsException", function(object) .print.OwsException(object))
setMethod("show", "OwsRange", function(object) .print.OwsRange(object))
setMethod("show", "SOS", function(object) .print.SOS(object))
setMethod("show", "SOS_1.0.0", function(object) .print.SOS_1.0.0(object))
setMethod("show", "SosFilter_Capabilities", function(object) .print.SosFilter_Capabilities(object))
setMethod("show", "SosObservationOffering", function(object) .print.SosObservationOffering(object))
setMethod("show", "SosContents", function(object) .print.SosContents(object))
setMethod("show", "SosEventTime", function(object) .print.SosEventTime(object))
setMethod("show", "SosFeatureOfInterest", function(object) .print.SosFeatureOfInterest(object))
setMethod("show", "SensorML", function(object) .print.SensorML(object))
setMethod("show", "SosGetObservation", function(object) .print.SosGetObservation(object))
setMethod("show", "SosGetObservationById", function(object) .print.SosGetObservationById(object))
setMethod("show", "SosDescribeSensor", function(object) .print.SosDescribeSensor(object))
setMethod("show", "SaSamplingPoint", function(object) .print.SaSamplingPoint(object))
setMethod("show", "SaSamplingSurface", function(object) .print.SaSamplingSurface(object))
setMethod("show", "SwePhenomenon", function(object) .print.SwePhenomenon(object))
setMethod("show", "SwePhenomenonProperty", function(object) .print.SwePhenomenonProperty(object))
setMethod("show", "SweCompositePhenomenon", function(object) .print.SweCompositePhenomenon(object))
setMethod("show", "SweTextBlock", function(object) .print.SweTextBlock(object))
setMethod("show", "OmObservationCollection", function(object) .print.OmObservationCollection(object))
setMethod("show", "OmObservation", function(object) .print.OmObservation(object))
setMethod("show", "OmObservationProperty", function(object) .print.OmObservationProperty(object))
setMethod("show", "GmlMeasure", function(object) .print.GmlMeasure(object))
setMethod("show", "OmMeasurement", function(object) .print.OmMeasurement(object))
setMethod("show", "GmlTimePosition", function(object) .print.GmlTimePosition(object))
setMethod("show", "GmlTimeInstant", function(object) .print.GmlTimeInstant(object))
setMethod("show", "GmlTimeInterval", function(object) .print.GmlTimeInterval(object))
setMethod("show", "GmlTimePeriod", function(object) .print.GmlTimePeriod(object))
setMethod("show", "GmlFeatureProperty", function(object) .print.GmlFeatureProperty(object))
setMethod("show", "GmlFeatureCollection", function(object) .print.GmlFeatureCollection(object))
setMethod("show", "GmlDirectPosition", function(object) .print.GmlDirectPosition(object))
setMethod("show", "GmlPoint", function(object) .print.GmlPoint(object))
setMethod("show", "GmlPointProperty", function(object) .print.GmlPointProperty(object))
setMethod("show", "GmlGeometry", function(object) .print.GmlGeometry(object))
setMethod("show", "GmlEnvelope", function(object) .print.GmlEnvelope(object))
setMethod("show", "TM_After", function(object) .print.TM_After(object))
setMethod("show", "TM_Before", function(object) .print.TM_Before(object))
setMethod("show", "TM_During", function(object) .print.TM_During(object))
setMethod("show", "TM_Equals", function(object) .print.TM_Equals(object))
setMethod("show", "OgcBBOX", function(object) .print.OgcBBOX(object))
setMethod("show", "OgcContains", function(object) .print.OgcContains(object))
setMethod("show", "OgcIntersects", function(object) .print.OgcIntersects(object))
setMethod("show", "OgcOverlaps", function(object) .print.OgcOverlaps(object))

################################################################################
# SUMMARY FUNCTIONS
summary.SOS = function(object, ...) {
	obj = list()
	obj[["class"]] = class(object)
	obj[["version"]] = sosVersion(object)
	obj[["url"]] = sosUrl(object)
	obj[["method"]] = sosMethod(object)
	obj[["title"]] = sosTitle(object)
	obj[["abstract"]] = sosAbstract(object)
	
	if(!is.null(sosTime(object)))
		obj[["time"]] = summary(sosTime(object))
	else obj[["time"]] = NA_character_
	
	obj[["offeringCount"]] = length(sosOfferingIds(object))
	obj[["procedureCount"]] = length(unlist(sosProcedures(object)))
	obj[["observedPropCount"]] = length(unlist(sosObservedProperties(object)))

	class(obj) = "summary.SOS"
	obj
}
setMethod("summary", "SOS", summary.SOS)

print.summary.SOS = function(x, ...) {
	cat(paste("Object of class ", x[["class"]], "\n", sep = ""))
	cat("[[version:]]\t")
	print(x[["version"]])
	cat("[[url:]]\t")
	print(x[["url"]])
	cat("[[title:]]\t")
	print(x[["title"]])
	cat("[[method:]]\t")
	print(x[["method"]])
	cat("[[abstract:]]\t")
	print(x[["abstract"]])
	cat("[[time:]]\t")
	print(x[["time"]])
	cat("[[offerings:]]\t")
	print(x[["offeringCount"]])
	cat("[[procedures:]]\t")
	print(x[["procedureCount"]])
	cat("[[observed properties:]]\t")
	print(x[["observedPropCount"]])
	invisible(x)
}

summary.SosObservationOffering = function(object, ...) {
	obj = list()
	obj[["class"]] = class(object)
	obj[["id"]] = sosId(object)
	obj[["name"]] = sosName(object)
	obj[["time"]] = summary(sosTime(object))
	obj[["bbox"]] = toString(sosBoundedBy(object))
	obj[["foiCount"]] = length(sosFeaturesOfInterest(object))
	obj[["procedureCount"]] = length(unlist(sosProcedures(object)))
	obj[["observedPropCount"]] = length(unlist(sosObservedProperties(object)))
	class(obj) = "summary.SosObservationOffering"
	obj
}
setMethod("summary", "SosObservationOffering", summary.SosObservationOffering)

print.summary.SosObservationOffering = function(x, ...) {
	cat(paste("Object of class ", x[["class"]], "\n", sep = ""))
	cat("[[id:]]\t\t")
	print(x[["id"]])
	cat("[[name:]]\t")
	print(x[["name"]])
	cat("[[time:]]\t")
	print(x[["time"]])
	cat("[[bbox:]]\t")
	print(x[["bbox"]])
	cat("[[fois:]]\t")
	print(x[["foiCount"]])
	cat("[[procs:]]\t")
	print(x[["procedureCount"]])
	cat("[[obsProps:]]\t")
	print(x[["observedPropCount"]])
	invisible(x)
}

summary.OwsRange = function(object, ...) {
	obj = list()
	obj[["class"]] = class(object)
	obj[["range"]] = paste(object@minimumValue, "-->", object@maximumValue)
	class(obj) = "summary.OwsRange"
	obj
}
setMethod("summary", "OwsRange", summary.OwsRange)
print.summary.OwsRange = function(x, ...) {
	print(x[["range"]])
	invisible(x)
}

summary.GmlTimePeriod = function(object, ...) {
	if(!is.null(object@begin) && !is.null(object@end)) {
		.s <- paste(toString(object@begin), "-->", toString(object@end))
	}
	else {
		.s <- paste(toString(object@beginPosition@time), "-->",
				toString(object@endPosition@time))
	}
	
	obj = list()
	obj[["class"]] = class(object)
	obj[["duration"]] = object@duration
	obj[["interval"]] = object@timeInterval
	obj[["beginEnd"]] = .s
	class(obj) = "summary.GmlTimePeriod"
	obj
}
setMethod("summary", "GmlTimePeriod", summary.GmlTimePeriod)
print.summary.GmlTimePeriod = function(x, ...) {
	print(x[["beginEnd"]])
	invisible(x)
}

summary.OmObservation = function(object, ...) {
	obj = list()
	obj[["class"]] = class(object)
	obj[["samplingTime"]] = length(object@samplingTime)
	obj[["procedureCount"]] = length(object@procedure)
	obj[["obsPropCount"]] = length(object@observedProperty)
	obj[["featureCount"]] = length(object@featureOfInterest)
	obj[["result"]] = summary(object@result)

	class(obj) = "summary.OmObservation"
	obj
}
setMethod("summary", "OmObservation", summary.OmObservation)
print.summary.OmObservation = function(x, ...) {
	cat(paste("Object of class ", x[["class"]], "\n", sep = ""))
	cat("[[samplingTime:]]\t")
	print(x[["samplingTime"]])
	cat("[[procedures:]]\t\t")
	print(x[["procedureCount"]])
	cat("[[obs. props:]]\t\t")
	print(x[["obsPropCount"]])
	cat("[[features:]]\t\t")
	print(x[["featureCount"]])
	cat("[[result summary:]]\n")
	print(x[["result"]])
	invisible(x)
}

summary.OmObservationCollection = function(object, ...) {
	obj = list()
	obj[["class"]] = class(object)
	obj[["memberCount"]] = length(object@members)
	obj[["boundedBy"]] = toString(object@boundedBy)
	obj[["procedureCount"]] = length(unique(unlist(sosProcedures(object))))
	obj[["obsPropCount"]] = length(unique(unlist(
							sosObservedProperties(object))))
	obj[["featureCount"]] = length(unique(unlist(sosFeatureIds(object))))
	
	class(obj) = "summary.OmObservationCollection"
	obj
}
setMethod("summary", "OmObservationCollection", summary.OmObservationCollection)
print.summary.OmObservationCollection = function(x, ...) {
	cat(paste("Object of class ", x[["class"]], "\n", sep = ""))
	cat("[[members:]]\t\t")
	print(x[["memberCount"]])
	cat("[[bounded by:]]\t\t")
	print(x[["boundedBy"]])
	cat("[[procedures:]]\t\t")
	print(x[["procedureCount"]])
	cat("[[obs. props:]]\t\t")
	print(x[["obsPropCount"]])
	cat("[[features:]]\t\t")
	print(x[["featureCount"]])
	invisible(x)
}


################################################################################
# utils
.addTabIndent <- function(str) {
	.s <- gsub(pattern = "\t",  replacement = "\t\t", x = str)
	return(.s)
}

