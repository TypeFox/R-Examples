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
# construction functions
#
OwsGetCapabilities <- function(
		service,
		acceptVersions,
		sections = sosDefaultGetCapSections,
		acceptFormats = sosDefaultGetCapAcceptFormats,
		updateSequence = c(as.character(NA)),
		owsVersion = sosDefaultGetCapOwsVersion,
		acceptLanguages = c(NA)) {
	if(owsVersion == "1.1.0") {
		if(!any(sapply(acceptLanguages, "is.na"), na.rm = TRUE))
			warning("Parameter 'acceptLanguages' is lost because it is not included in 1.1.0!")
		new("OwsGetCapabilities_1.1.0",
				request = sosGetCapabilitiesName,
				version = "1.1.0",
				service = service,
				acceptVersions = acceptVersions, sections = sections,
				acceptFormats = acceptFormats, updateSequence = updateSequence)
	}
	else if(owsVersion == "2.0.0") {
		new("OwsGetCapabilities_2.0.0",
				request = sosGetCapabilitiesName,
				version = "2.0.0",
				service = service,
				acceptVersions = acceptVersions, sections = sections,
				acceptFormats = acceptFormats, updateSequence = updateSequence,
				acceptLanguages = acceptLanguages)
	}
	else {
		new("OwsGetCapabilities",
				request = sosGetCapabilitiesName, version = "NONE",
				service = service, acceptVersions = acceptVersions,
				owsVersion = owsVersion)
	}
}

OwsCapabilities <- function(
		version, 
		updateSequence = NA,
		owsVersion = sosDefaultGetCapOwsVersion,
		identification = NULL,
		provider = NULL,
		operations = NULL,
		contents = NULL,
		languages = NULL) {
	if(owsVersion == "1.1.0") {
		if(!is.na(languages))
			warning("Parameter 'languages' is lost because it is not included in 1.1.0!")
		new("OwsCapabilities_1.1.0",
			version = version, updateSequence = updateSequence,
			owsVersion = owsVersion, identification = identification,
			provider = provider, operations = operations,
			contents = contents)
	}
	else if(owsVersion == "2.0.0") {
		new("OwsCapabilities_2.0.0",
			version = version, updateSequence = updateSequence,
			owsVersion = owsVersion, identification = identification,
			provider = provider, operations = operations,
			contents = contents, languages = languages)
	}
	else {
		new("OwsCapabilities",
			version = version, updateSequence = updateSequence,
			owsVersion = owsVersion)
	}	
}

OwsServiceIdentification <- function(serviceType, serviceTypeVersion,
		profile = c(NA), title, abstract = c(NA), keywords = c(NA),
		fees = as.character(NA), accessConstraints = c(NA)) {
	new("OwsServiceIdentification",
			serviceType = serviceType, serviceTypeVersion = serviceTypeVersion,
			profile = profile, title = title, abstract = abstract,
			keywords = keywords, fees = fees,
			accessConstraints = accessConstraints)
}

OwsServiceProvider <- function(providerName, providerSite = as.character(NA),
		serviceContact = xmlNode(NA)) {
	new("OwsServiceProvider", providerName = providerName,
			providerSite = providerSite, serviceContact = serviceContact)
}

OwsOperationsMetadata <- function(operations, parameters = list(NA),
	constraints = list(NA), extendedCapabilities = xmlNode(NA)) {
	new("OwsOperationsMetadata", operations = operations,
		parameters = parameters, constraints = constraints,
		extendedCapabilities = extendedCapabilities)
}

OwsOperation <- function(name, DCPs, parameters = list(NA),
		constraints = list(NA), metadata = list(NA)) {
	new("OwsOperation", name = name, DCPs = DCPs, parameters = parameters,
		constraints = constraints, metadata = metadata)
}

OwsContents <- function(xmlNode) {
	new("OwsContents", xml = xmlNode)
}

OwsExceptionReport <- function(version, lang = as.character(NA),
		exceptions = list(NA)) {
	new("OwsExceptionReport", version = version, lang = lang,
			exceptions = exceptions)
}

OwsException <- function(exceptionCode, exceptionText = c(),
		locator = as.character(NA)) {
	new("OwsException", exceptionCode = exceptionCode,
			exceptionText = exceptionText, 
			locator = locator)
}

OwsRange <- function(minimumValue = as.character(NA),
		maximumValue = as.character(NA), rangeClosure = as.character(NA),
		spacing = as.character(NA)) {
	new("OwsRange", minimumValue = minimumValue, maximumValue = maximumValue,
			rangeClosure = rangeClosure, spacing = spacing)
}



################################################################################
# checking of operations before they are sent out
#
setMethod(f = "checkRequest",
		signature = signature(service = "SOS",
				operation = "OwsGetCapabilities_1.1.0",
				verbose = "logical"),
		def = function(service, operation, verbose) {
			
			# TODO implement checkRequest for OwsGetCapabilities
			
			return(TRUE)
		})
setMethod(f = "checkRequest",
		signature = signature(service = "SOS",
				operation = "OwsGetCapabilities_2.0.0",
				verbose = "logical"),
		def = function(service, operation, verbose) {
			
			# TODO implement checkRequest for OwsGetCapabilities
			
			return(TRUE)
		})


################################################################################
# helper methods
#

# 
# to add (possible) multiple values in kvp
#
.kvpKeyAndValues <- function(key, values) {
	if(is(values, "vector")) {
		.values <- sapply(values, .kvpEscapeSpecialCharacters)
		valueList <- paste(.values, collapse = ",")
		return(paste(key, valueList, sep = "="))
	}
	else {
		return(paste(key, .kvpEscapeSpecialCharacters(x = values), sep = "="))
	}
}

#
# Method to excape characters within values (!) of a parameter. This function
# cannot be called on the whole request string!
#
# See http://www.ietf.org/rfc/rfc2396.txt and 
# http://www.oostethys.org/best-practices/best-practices-get and
# http://www.opengeospatial.org/standards/common (Section 11.3)
# and maybe also http://www.blooberry.com/indexdot/html/topics/urlencoding.htm
#
# Special character  	Escaped encoding
# :					 	%3A
# / 					%2F
# # 					%23
# ? 					%3F
# = 					%3D
#
.kvpEscapeSpecialCharacters <- function(x) {
	.escaped <- gsub(x = x, pattern = ":", replacement = "%3A")
	#.escaped <- gsub(.escaped, pattern = "/", replacement = "%2F")
	.escaped <- gsub(x = .escaped, pattern = "#", replacement = "%23")
	.escaped <- gsub(x = .escaped, pattern = "\\?", replacement = "%3F")
	.escaped <- gsub(x = .escaped, pattern = "=", replacement = "%3D")
	.escaped <- gsub(x = .escaped, pattern = "&", replacement = "%26")
	.escaped <- gsub(x = .escaped, pattern = ",", replacement = "%2C")
	.escaped <- gsub(x = .escaped, pattern = "\\+", replacement = "%2B")
	.escaped <- gsub(x = .escaped, pattern = "@", replacement = "%40")
	return(.escaped)
}

################################################################################
# kvp encoding
#
setMethod(f = "encodeRequestKVP", "OwsGetCapabilities",
		def = function(obj, sos, verbose = FALSE) {
			.sosEncodeRequestKVPGetCapabilities(obj, verbose)
		})
.sosEncodeRequestKVPGetCapabilities <- function(obj, verbose = FALSE) {
	.service <- paste(
			"service",
			.kvpEscapeSpecialCharacters(x = obj@service),
			sep = "=")
	.request <- paste(
			"request",
			.kvpEscapeSpecialCharacters(x = obj@request),
			sep = "=")
	
	.kvpString <- paste(.service, .request, sep = "&")
	
	return(.kvpString)
}

setMethod(f = "encodeRequestKVP", "OwsGetCapabilities_1.1.0",
		function(obj, sos, verbose = FALSE) {
			.sosEncodeRequestKVPGetCapabilities_1.1.0(obj, verbose)
		})
.sosEncodeRequestKVPGetCapabilities_1.1.0 <- function(obj, verbose = FALSE) {
	.mandatory <- .sosEncodeRequestKVPGetCapabilities(obj, verbose)
	
	.optionals = ""
	if( !is.na(obj@acceptVersions)) {
		.optionals <- paste(.optionals, .kvpKeyAndValues("acceptVersions",
						obj@acceptVersions), sep = "&")
	}
	
	if(!any(sapply(obj@sections, "is.na"), na.rm = TRUE)) {
		.optionals <- paste(.optionals, .kvpKeyAndValues("sections",
						obj@sections), sep = "&")
	}
	
	if( !is.na(obj@updateSequence)) {
		.optionals <- paste(.optionals, .kvpKeyAndValues("updateSequence",
						obj@updateSequence), sep = "&")
	}
	
	if(!any(sapply(obj@acceptFormats, "is.na"), na.rm = TRUE)) {
		.optionals <- paste(.optionals, .kvpKeyAndValues("acceptFormats",
						obj@acceptFormats), sep = "&")
	}
	
	.kvpString <- paste(.mandatory, .optionals, sep = "")
	
	if(verbose) cat("[.sosEncodeRequestKVPGetCapabilities_1.1.0] done: ",
				.kvpString, "\n")
	
	return(.kvpString)
}

setMethod(f = "encodeRequestKVP", "OwsGetCapabilities_2.0.0",
		function(obj, sos, verbose = FALSE) {
			.sosEncodeRequestKVPGetCapabilities_2.0.0(obj, verbose)
		})
.sosEncodeRequestKVPGetCapabilities_2.0.0 <- function(obj, verbose = FALSE) {
	.kvpString <- .sosEncodeRequestKVPGetCapabilities_1.1.0(obj)
	
	if(!any(sapply(obj@acceptLanguages, "is.na"), na.rm = TRUE)) {
		.kvpString <- paste(.kvpString, .kvpKeyAndValues("acceptLanguages",
						obj@acceptLanguages), sep = "&")
	}
	
	if(verbose) cat("[.sosEncodeRequestKVPGetCapabilities_2.0.0] done: ",
				.kvpString, "\n")
	
	return(.kvpString)
}

################################################################################
# XML encoding
#
setMethod("encodeRequestXML", "OwsGetCapabilities_1.1.0", 
		function(obj, sos, verbose = FALSE) {
			if(verbose) {
				cat("[encodeRequestXML]", class(obj), "\n")
			}
			
			return(.sosEncodeRequestXMLOwsGetCapabilities_1.1.0(obj))
		}
)
.sosEncodeRequestXMLOwsGetCapabilities_1.1.0 <- function(obj) {
	.xmlDoc <- xmlNode(name = sosGetCapabilitiesName,
			namespace = sosNamespacePrefix,
			namespaceDefinitions = c(.sosNamespaceDefinitionsForAll,
					.sosNamespaceDefinitionsGetCap),
			attrs=c(.xsiSchemaLocationAttribute,
					service=obj@service))
	
	# optional:
	if( !is.na(obj@acceptVersions)) {
		.acceptVersions <- xmlNode(name = "AcceptVersions",
				namespace = owsNamespacePrefix)
		.acceptVersions$children <- lapply(
				obj@acceptVersions, "xmlNode", name="ows:Version")
		.xmlDoc <- addChildren(node = .xmlDoc, kids = list(.acceptVersions))
	}
	
	if(!any(sapply(obj@sections, "is.na"), na.rm = TRUE)) {
		.sections <- xmlNode("ows:Sections")
		.sections$children <- lapply(obj@sections, "xmlNode", name="Section",
				namespace = owsNamespacePrefix)
		.xmlDoc <- addChildren(node = .xmlDoc, kids = list(.sections))
	}
	
	if( !is.na(obj@updateSequence)) {
		.xmlDoc <- addAttributes(.xmlDoc, updateSequence = obj@updateSequence)
	}
	
	if(!any(sapply(obj@acceptFormats, "is.na"), na.rm = TRUE)) {
		.acceptFormats <- xmlNode(name = "AcceptFormats",
				namespace = owsNamespacePrefix)
		.acceptFormats$children <- lapply(
				obj@acceptFormats, "xmlNode", name="ows:OutputFormat")
		.xmlDoc <- addChildren(node = .xmlDoc, kids = list(.acceptFormats))
	}
	
	return(.xmlDoc)
}

#
#
#
setMethod("encodeRequestXML", "OwsGetCapabilities_2.0.0", 
		function(obj, sos, verbose = FALSE) {
			if(verbose) {
				cat("[encodeRequestXML]", class(obj), "\n")
			}
			
			return(.sosEncodeRequestXMLOwsGetCapabilities_2.0.0(obj))
		}
)
.sosEncodeRequestXMLOwsGetCapabilities_2.0.0 <- function(obj) {
	.xmlDoc <- xmlNode(name = sosGetCapabilitiesName,
			namespace = sosNamespacePrefix,
			namespaceDefinitions = c(.sosNamespaceDefinitionsForAll,
					.sosNamespaceDefinitionsGetCap),
			attrs=c(.xsiSchemaLocationAttribute,
					service=obj@service))
	
	# optional:
	if( !is.na(obj@acceptVersions)) {
		.acceptVersions <- xmlNode(name = "AcceptVersions",
				namespace = owsNamespacePrefix)
		.acceptVersions$children <- lapply(
				obj@acceptVersions, "xmlNode", name = "ows:Version")
		.xmlDoc <- addChildren(node = .xmlDoc, kids = list(.acceptVersions))
	}
	
	if(!any(sapply(obj@sections, "is.na"), na.rm = TRUE)) {
		.sections <- xmlNode("ows:Sections")
		.sections$children <- lapply(obj@sections, "xmlNode", name = "Section",
				namespace = owsNamespacePrefix)
		.xmlDoc <- addChildren(node = .xmlDoc, kids = list(.sections))
	}
	
	if( !is.na(obj@updateSequence)) {
		.xmlDoc <- addAttributes(.xmlDoc, updateSequence = obj@updateSequence)
	}
	
	if(!any(sapply(obj@acceptFormats, "is.na"), na.rm = TRUE)) {
		.acceptFormats <- xmlNode(name = "AcceptFormats",
				namespace = owsNamespacePrefix)
		.acceptFormats$children <- lapply(
				obj@acceptFormats, "xmlNode", name="ows:OutputFormat")
		.xmlDoc <- addChildren(node = .xmlDoc, kids = list(.acceptFormats))
	}
	
	if(!any(sapply(obj@acceptLanguages, "is.na"), na.rm = TRUE)) {
		.acceptLanguages <- xmlNode(name = "AcceptLanguages",
				namespace = owsNamespacePrefix)
		.acceptLanguages$children <- lapply(
				obj@acceptLanguages, "xmlNode", name="ows:Language")
		.xmlDoc <- addChildren(node = .xmlDoc, kids = list(.acceptLanguages))
	}
	
	return(.xmlDoc)
}


################################################################################
# SOAP encoding
setMethod("encodeRequestSOAP", "OwsGetCapabilities", 
		function(obj, sos, verbose = FALSE) {
			if(verbose) {
				cat("ENCODE SOAP ", class(obj), "\n")
			}
			
			stop("Function not implemented yet...")
		}
)
