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

#
# virtual parent class for all service requests etc., loosely based on
# commonalities of OWS 1.1.0 and OWS 2.0.0
#
setClass("OwsServiceOperation",
		representation(
				service = "character",	# service type identifier, e.g. "WMS"
				request = "character",	# operation name, e.g. "GetMap"
				version = "character"),	# version of the operation
		prototype = list(service = as.character(NA), request= as.character(NA),
				version = as.character(NA)),
		contains = "VIRTUAL",
		validity = function(object) {
			#print("Entering validation: OwsServiceOperation")
			# TODO implement validity function
			
			if(is.na(object@service))
				return("service parameter must be given")
			if(is.na(object@request))
				return("request parameter must be given")
			if(is.na(object@version))
				return("version parameter must be given")
			
			return(TRUE)
		}
)

#
# Mandatory parameters, see OWS Common 2.0, OGC 06-121r3
#
setClass("OwsGetCapabilities",
		representation(acceptVersions = "character", owsVersion = "character"),
		prototype = list(service = as.character(NA),
				request = sosGetCapabilitiesName,
				acceptVersions = as.character(NA),
				owsVersion = as.character(NA)),
		contains = "OwsServiceOperation",
		validity = function(object) {
			#print("Entering validation: OwsGetCapabilities")
			# TODO implement validity function
			
			# service and request need to be there
			if(is.na(object@service))
				return("service parameter must be given")
			if(is.na(object@acceptVersions))
				return("acceptVersions vector must be given")
			# if acceptVersion is there, it hast to be in a certain format, see ows common
			
			# owsVersion has to be one of 1.1.0 or 2.0.0
			.allowedOwsVersions <- c("1.1.0", "2.0.0")
			if(!any(sapply(.allowedOwsVersions, "==", object@owsVersion)))
				return(paste("owsVersion must be one of", paste(.allowedOwsVersions, collapse = ", ")))
			
			return(TRUE)
		}
)

#
# See OWS Common 1.1.0, OGC 06-121r3
#
setClass("OwsGetCapabilities_1.1.0",
		representation(sections = "vector", acceptFormats = "vector",
				updateSequence = "vector"),
		prototype = list(owsVersion = "1.1.0"),
		contains = "OwsGetCapabilities",
		validity = function(object) {
			#print("Entering validation: OwsGetCapabilities_1.1.0")
			# TODO implement validity function
			
			# service and request need to be there
			if(is.na(object@service))
				return("service parameter must be given")
			if(is.na(object@acceptVersions))
				return("acceptVersions vector must be given")
			
			# if acceptVersion is there, it hast to be in a certain format, see ows common
			
			return(TRUE)
		}
)

#
# See OWS Common 2.0, OGC 06-121r9
#
setClass("OwsGetCapabilities_2.0.0",
		representation(acceptLanguages = "vector"),
		prototype = list(acceptLanguages = c(NA),
				owsVersion = "2.0.0"),
		contains = "OwsGetCapabilities_1.1.0",
		validity = function(object) {
			#print("Entering validation: OwsGetCapabilities_2.0.0")
			# TODO implement validity function
			
			# service and request need to be there
			if(is.na(object@service))
				return("service parameter must be given")
			if(is.na(object@acceptVersions))
				return("acceptVersions vector must be given")
			
			# TODO add check: if acceptVersion is there, it hast to be in a certainformat, see OWS Common...
			
			return(TRUE)
		}
)

################################################################################
# Capabilities Content:

#
# See OGC 06-121r3, clause 7.4.6 
#
setClass("OwsOperationsMetadata",
		representation(operations = "list", parameters = "list",
				constraints = "list",
				extendedCapabilities = "XMLAbstractNode"),
		prototype = list(operations = list(NA)),
		validity = function(object) {
			#print("Entering validation: OwsOperationsMetadata")
			# TODO implement validity function
			
			# operations must all be of class OwsOperation and at least one must
			# be there
			
			return(TRUE)
		}
)
setClassUnion(name = "OwsOperationsMetadataOrNULL",
		members = c("OwsOperationsMetadata", "NULL"))

#
# See OGC 06-121r3, clause 7.4.6
#
setClass("OwsOperation",
		representation(name = "character", DCPs = "list",
			parameters = "list", constraints = "list",
			metadata = "list"),
		prototype = list(name = as.character(NA), DCPs = list(NA)),
		validity = function(object) {
			#print("Entering validation: OwsOperation")
			# TODO implement validity function
			
			# name is mandatory
			# one dcp is mandatory (at the moment only http anyway, but more possible)
			# other parameters are optional
			
			return(TRUE)
		}
)

#
# See OGC 06-121r3, clause 7.4.4 or OGC 06-121r9, clause 7.4.4 (only changes is
# the data type of the parameter 'serviceType' that changed to URN instead of
# character string type, and that does not require special handling here).
#
setClass("OwsServiceIdentification",
		representation(serviceType = "character",
				serviceTypeVersion = "vector",
				profile = "vector",
				title = "vector",
				abstract = "vector",
				keywords = "vector",
				fees = "character",
				accessConstraints = "vector"),
		prototype = list(),
		validity = function(object) {
			#print("Entering validation: OwsServiceIdentification")
			# TODO implement validity function
			
			# mandatory elements: serviceType, serviceTypeVersion, title
			
			return(TRUE)
		}
)
setClassUnion(name = "OwsServiceIdentificationOrNULL",
		members = c("OwsServiceIdentification", "NULL"))

#
# See OGC 06-121r3, clause 7.4.5 or OGC 06-121r9, clause 7.4.4 (no changes
# between the two versions).
#
setClass("OwsServiceProvider",
		representation(providerName = "character", providerSite = "character",
				serviceContact = "XMLAbstractNode"),
		prototype = list(providerName = as.character(NA)),
		validity = function(object) {
			#print("Entering validation: OwsServiceProvider")
			
			#TODO implement validity function, providerName is mandatory
			
			return(TRUE)
		}
)
setClassUnion(name = "OwsServiceProviderOrNULL",
		members = c("OwsServiceProvider", "NULL"))

#
# See OGC 06-121r3, clause 7.4.8
#
setClass("OwsContents",
		representation(xml = "XMLAbstractNode"),
		prototype = list(xml = NULL),
		validity = function(object) {
			#print("Entering validation: OwsContents")
			return(TRUE)
		}
)
setClassUnion(name = "OwsContentsOrNULL",
		members = c("OwsContents", "NULL"))

################################################################################
# Capabilities documents:

#
# Mandatory parameters, see OWS Common 2.0, OGC 06-121r3
#
setClass("OwsCapabilities",
		representation(version = "character", 
			updateSequence = "character",
			owsVersion = "character"),
		prototype = list(version = as.character(NA),
			updateSequence = as.character(NA),
			owsVersion = as.character(NA)),
		validity = function(object) {
			#print("Entering validation: OwsCapabilities")
			# TODO implement validity function
			
			# version is mandatory
			if(is.na(object@version))
				return("version parameter must be given")
			
			# owsVersion has to be one of 1.1.0 or 2.0.0
			.allowedOwsVersions <- c("1.1.0", "2.0.0")
			if(!any(sapply(.allowedOwsVersions, "==", object@owsVersion)))
				return(paste("owsVersion must be one of", paste(.allowedOwsVersions, collapse = ", ")))
			
			return(TRUE)
		}
)

#
# See OWS Common 1.1.0, OGC 06-121r3
#
setClass("OwsCapabilities_1.1.0",
		representation(identification = "OwsServiceIdentificationOrNULL",
				provider = "OwsServiceProviderOrNULL",
				operations = "OwsOperationsMetadataOrNULL",
				contents = "OwsContentsOrNULL"),
		prototype = list(owsVersion = "1.1.0",
			identification = NULL, provider = NULL,
			operations = NULL, contents = NULL),
		contains = "OwsCapabilities",
		validity = function(object) {
			#print("Entering validation: OwsCapabilities_1.1.0")
			# TODO implement validity function
			return(TRUE)
		}
)

#
# See OWS Common 2.0, OGC 06-121r9
# languages elements are represented by character only!
#
setClass("OwsCapabilities_2.0.0",
		representation(languages = "XMLAbstractNode"),
		prototype = list("GetCapabilities",
				languages = xmlNode(NA),
				owsVersion = "2.0.0"),
		contains = "OwsCapabilities_1.1.0",
		validity = function(object) {
			#print("Entering validation: OwsCapabilities_2.0.0")
			# TODO implement validity function
			return(TRUE)
		}
)


################################################################################
# Service Exceptions:

#
# See OWS Common 1.0 and 2.0, OGC 06-121r3 respectively r9
#
# (HTTP Status Codes (section 8.6) and SOAP encoding (8.7) as described in
# version 2.0.0 are not included here.)
#
setClass("OwsExceptionReport",
		representation(version = "character", lang = "character",
				exceptions = "list"),
		prototype = list(version = as.character(NA), lang = as.character(NA),
				exceptions = list()),
		validity = function(object) {
			#print("Entering validation: OwsExceptionReport")
			# TODO implement validity function
			
			# version needs to be there
			if(is.na(object@version))
				return("version parameter must be given")
			
			return(TRUE)
		}
)

#
#
#
setClass("OwsException",
		representation(exceptionCode = "character", exceptionText = "vector",
				locator = "character"),
		prototype = list(code = as.character(NA)),
		validity = function(object) {
			#print("Entering validation: OwsException")
			# TODO implement validity function
			
			# version needs to be there
			if(is.na(object@exceptionCode))
				return("exceptionCode parameter must be given")
			
			return(TRUE)
		}
)

#
#
#
setClass("OwsRange",
		representation(minimumValue = "character", maximumValue = "character",
				rangeClosure = "character", spacing = "character"),
		prototype = list(),
		validity = function(object) {
			#print("Entering validation: OwsRange")
			# TODO implement validity function
			
			# all elements are optional!
			
			# closure: one of closed, open, open-closed, closed-open, see owsDomainType.xsd
			return(TRUE)
		}
)

