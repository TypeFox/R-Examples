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
#
#
parseOwsOperation <- function(obj) {
	.name <- xmlGetAttr(obj, "name")
	
	.dcpsXML <- .filterXmlChildren(obj, owsDCPName)
	.dcps <- list()
	for(.dcp in .dcpsXML) {
		.http <- .dcp[[owsHTTPName]]
		.endpoints <- c(
				.filterXmlChildren(.http, owsGetName),
				.filterXmlChildren(.http, owsPostName))
		
		for(.ep in .endpoints) {
			.newEndpoint <- list(xmlGetAttr(.ep, "href"))
			names(.newEndpoint) <- xmlName(.ep)
			.dcps <- c(.dcps, .newEndpoint)
		}
	}
	
	.parametersXML <- .filterXmlChildren(obj, owsParameterName)
	.parameters = list()
	.names = list()
	
	if(length(.parametersXML) > 0) {
		for(.p in .parametersXML) {
			.allowedValues <- NULL
			.ranges <- NULL
			.allowedValuesAndRanges <- NULL
			
			# check for ows:AnyValue	
			if(length(.p[owsAnyValueName]) > 0)
				.allowedValuesAndRanges = list(owsAnyValueName)
			else {
				# list of allowed values
				.xpathAllowedValues <- paste("./", owsNamespacePrefix, ":",
						owsAllowedValuesName, "/", owsNamespacePrefix, ":",
						owsValueName, sep = "")
				.allowedValues <- lapply(
						getNodeSet(doc = .p, path = .xpathAllowedValues,
								namespaces = .owsNamespace),
						xmlValue)
				# list of ranges
				.xpathRanges <- paste("./", owsNamespacePrefix, ":",
						owsAllowedValuesName, "/", owsNamespacePrefix, ":",
						owsRangeName, sep = "")
				.ranges <-  sapply(
						getNodeSet(
								.p,
								.xpathRanges,
								.owsNamespace
						),
						parseOwsRange)
				.allowedValuesAndRanges <- c(.allowedValues, .ranges)
			}
			
#			cat("[", .name, "] Adding to parameters list for",
#					xmlGetAttr(.p, "name"), ":",
#					toString(.allowedValuesAndRanges), "\n")
			
			.names <- c(.names, xmlGetAttr(.p, "name"))
			.parameters[[length(.parameters) + 1]] <- .allowedValuesAndRanges
			# the following does NOT work as it recursively concatenates the
			# lists: .parameters <- c(.parameters, .allowedValuesAndRanges)
		}
		
		names(.parameters) <- .names
	}
	
	if(any(sapply(names(obj), "==", owsConstraintName)))
		warning("constraint elements are NOT processed!")
	.constraints = list(NA)
	
	if(any(sapply(names(obj), "==", owsMetadataName)))
		warning("metadata elements are NOT processed!")
	.metadata = list(NA)
	
	.op <- OwsOperation(name = .name, DCPs = .dcps,
			parameters = .parameters, constraints = .constraints,
			metadata = .metadata)
	return(.op)
}

#
# method for parsing an ows:ExceptionReport.
#
parseOwsExceptionReport <- function(obj, verbose = FALSE) {
	if(verbose) cat("[parseOwsExceptionReport] Starting ...")
	.docRoot <- xmlRoot(obj)
	## print(.docRoot)

	.version <- xmlGetAttr(node = .docRoot, name = "version")
	.lang <- xmlGetAttr(node = .docRoot, name = "lang", default = NA_character_)
	
	# remove all elements from docRoot that are not 'Exception'
	# could probably be done nicer with subsetting, but indexing with wildcards or similar (... xmlChildren()[[]] ...) did not work.
	.children <- xmlChildren(.docRoot) 
	.exceptionsXML <- list()
	for (x in .children) {
		if(xmlName(x) == owsExceptionName)
			.exceptionsXML = c(.exceptionsXML, x)
		# else print(xmlName(x))
	}
	
	.exceptions = sapply(.exceptionsXML, parseOwsException)
	if(verbose) cat("[parseOwsExceptionReport]", length(.exceptions),
				"exceptions.")
	
	.report <- OwsExceptionReport(version = .version, lang = .lang, exceptions = .exceptions)
	
	return(.report)
}

#
# parsing a single xml node that is an ows:Exception
#
parseOwsException <- function(obj) {
	#	print("parsing e!")
	.code <- xmlGetAttr(node = obj, name = "exceptionCode")
	.locator <- xmlGetAttr(node = obj, name = "locator",
			default = NA_character_)
	
	if(!is.na(xmlChildren(obj)[owsExceptionTextName]))
		.text <- xmlValue(xmlChildren(obj)[[owsExceptionTextName]])
	else .text <- as.character(NA)
	
	.exception <- OwsException(exceptionCode = .code, 
			exceptionText = .text,
			locator = .locator)
	
	return(.exception)
}

#
#
#
parseOwsServiceIdentification <- function(obj) {
#	print("parsing ows service identification!")
	
	.children <- xmlChildren(obj)
	.serviceType <- sapply(.filterXmlChildren(obj, owsServiceTypeName),
			xmlValue)
	.serviceTypeVersion <- sapply(.filterXmlChildren(obj,
					owsServiceTypeVersionName),
			xmlValue)
	.title <- sapply(.filterXmlChildren(obj, owsTitleName),
			xmlValue)
	
	# optional:
	if(!is.na(xmlChildren(obj)[owsProfileName]))
		.profile <- lapply(.filterXmlChildren(obj, owsProfileName), xmlValue)
	else .profile <- c(NA)
	
	if(!is.na(xmlChildren(obj)[owsAbstractName]))
		.abstract <- lapply(.filterXmlChildren(obj, owsAbstractName), xmlValue)
	else .abstract <- c(NA)
	
	if(!is.na(xmlChildren(obj)[owsKeywordsName])) {
		.keywordLists <- .filterXmlChildren(obj, owsKeywordsName)
		.keywords <- c(lapply(.keywordLists, FUN = xmlToList), recursive = TRUE)
		.keywords <- lapply(.keywords, gsub, pattern = "^[[:space:]]+|[[:space:]]+$",
				replacement = "") # http://finzi.psych.upenn.edu/R/Rhelp02a/archive/40714.html
	}
	else .keywords <- c(NA)
	
	if(!is.na(xmlChildren(obj)[owsFeesName]))
		.fees <- paste(sapply(.filterXmlChildren(obj, owsFeesName), xmlValue))
	else .fees <- as.character(NA)
	
	if(!is.na(xmlChildren(obj)[owsAccessConstraintsName]))
		.accessConstraints <- lapply(.filterXmlChildren(obj,
						owsAccessConstraintsName),
				xmlValue)
	else .accessConstraints <- c(NA)
	
	.si <- OwsServiceIdentification(serviceType =  .serviceType,
			serviceTypeVersion = .serviceTypeVersion, profile = .profile,
			title = .title, abstract = .abstract, keywords = .keywords,
			fees = .fees, accessConstraints = .accessConstraints)
}

#
#
#
parseOwsServiceProvider <- function(obj) {
	#print("parsing ows service provider!")
	.name <- xmlValue(obj[[owsProviderNameName]])
	
	# optional:
	if(!is.null(xmlChildren(obj)[[owsProviderSiteName]]))
		.site <- xmlGetAttr(node = obj[[owsProviderSiteName]],
				name = "href", default = as.character(NA))
	else .site <- as.character(NA)
	
	if(!is.null(xmlChildren(obj)[[owsServiceContactName]])) {
		.contact <- obj[[owsServiceContactName]]
		.sp <- OwsServiceProvider(providerName = .name, providerSite = .site,
				serviceContact = .contact)
	}
	else .sp <- OwsServiceProvider(providerName = .name, providerSite = .site)
	
	return(.sp)
}

#
# all elements are optional
#
parseOwsRange <- function(obj) {
	.children <- xmlChildren(obj)
	
	.minimumXml <- .children[[owsMinimumValueName]]
	if(is.null(.minimumXml)) {
		.minimum <- as.character(NA)
	} else {
		.minimum <- xmlValue(.minimumXml)
	}
	
	
	.maximumXml <- .children[[owsMaximumValueName]]
	if(is.null(.maximumXml)) {
		.maximum <- as.character(NA)
	} else {
		.maximum <- xmlValue(.maximumXml)
	}
	
	.closure <- xmlGetAttr(node = obj, name = "rangeClosure")
	if(is.null(.closure)) {
		.closure <- as.character(NA)
	}
	
	.spacingXml <- .children[[owsSpacingName]]
	if(is.null(.spacingXml)) {
		.spacing <- as.character(NA)
	} else {
		.spacing <- xmlValue(.spacingXml)
	}
	
	.range <- OwsRange(minimumValue = .minimum, maximumValue = .maximum,
			rangeClosure = .closure, spacing = .spacing)
	
	return(.range)
}


#
# If includeNamed is TRUE, this function returns a list of all children with the
# given name, of the given element. If includeNamed is FALSE, it contains all
# children of the given node not matching the given name.
#
.filterXmlChildren <- function(node, childrenName, includeNamed = TRUE,
		verbose = FALSE) {
	.temp <- xmlChildren(node)
	
	if(verbose) {
		cat("[.filterXmlChildren] Children:\n")
		print(.temp)
	}
	
	.filtered <- c()
	.names <- c()
	for (.x in .temp) {
		if(xmlName(.x) == childrenName && includeNamed) {
			.filtered <- c(.filtered, .x)
			if(verbose) cat("[.filterXmlChildren] Added", xmlName(.x), "\n")
			.names <- c(.names, xmlName(.x))
		}
		else if(!includeNamed && xmlName(.x) != childrenName) {
			.filtered <- c(.filtered, .x)
			if(verbose) cat("[.filterXmlChildren] Added", xmlName(.x), "\n")
			.names <- c(.names, xmlName(.x))
		}
	}
	names(.filtered) <- .names
	rm(.temp)
	rm(.names)
	return(.filtered)
}

.filterXmlOnlyNoneTexts <- function(node) {
	.filterXmlChildren(
			node = node,
			childrenName = xmlTextNodeName, includeNamed = FALSE)
}
