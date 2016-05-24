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
# Created: 2010-09-15                                                          #
# Project: sos4R - visit the project web page, http://www.nordholmen.net/sos4r #
#                                                                              #
################################################################################

#
# optimized for 52N SOS, that means only options used there in OMEncoder are
# handled here.
#
# For example, swe:elementCount can also have an attribut ref, but this is not
# checked here, and swe:Count is actually a swe:AbstractDataComponentType, but
# here it is just looked for a child element swe:value.
#
parseDataArray <- function(obj, sos, verbose = FALSE) {
	.elementCount <-  xmlValue(obj[["elementCount"]][["Count"]][["value"]])
	if(verbose) cat("[parseDataArray] Parsing DataArray with", .elementCount,
				"elements.\n")
	
	.eTParser <- sosParsers(sos)[[sweElementTypeName]]
	.fields <- .eTParser(obj = obj[[sweElementTypeName]], verbose = verbose)
	
	if(verbose) cat("[parseDataArray]  Parsed field descriptions:",
				toString(.fields), "\n")
	
	.encParser <- sosParsers(sos)[[sweEncodingName]]
	.encoding <- .encParser(obj = obj[[sweEncodingName]],
		verbose = verbose)
	
	if(verbose) cat("[parseDataArray]  Parsed encoding description:",
				toString(.encoding), "\n")
	
	.valParser <- sosParsers(sos)[[sweValuesName]]
	.values <- .valParser(values = obj[[sweValuesName]], fields = .fields,
			encoding = .encoding, sos = sos, verbose = verbose)
	
	return(.values)
}


#
# values is XML and encoding holds a SweTextBlock with the required separators.
#
parseValues <- function(values, fields, encoding, sos, verbose = FALSE) {
	if(verbose) cat("[parseValues] Parsing swe:values using", toString(encoding), "and",
				length(fields), "fields:", toString(names(fields)), "\n")
	if(!inherits(encoding, "SweTextBlock")) {
		stop("Handling for given encoding not implemented!")
	}
	
	.converters <- sosDataFieldConverters(sos)
	
	.blockLines <- strsplit(x = xmlValue(values),
			split = encoding@blockSeparator)
	.tokenLines <- sapply(.blockLines, strsplit,
			split = encoding@tokenSeparator)
	
	if(verbose)
		cat("[parseValues] Parsing values from lines: ", toString(.tokenLines), "\n")
	
	# data frame of correct length to be able to use cbind for first column
	.tempId = "tempID"
	.data <- data.frame(seq(1,length(.tokenLines)))
	names(.data) <- .tempId
	
	# do following for all fields
	.fieldCount <- length(fields)
	for (.currentFieldIdx in seq(1,.fieldCount)) {
		if(verbose)
			cat("[parseValues] Processing field index", .currentFieldIdx , "of", .fieldCount,"\n")
		
		# create list for each variable
		.currentValues <- sapply(.tokenLines, "[[", .currentFieldIdx)
		if(verbose)
			cat("[parseValues] Current values: ", toString(.currentValues), "\n")
		.currentField <- fields[[.currentFieldIdx]]
		
		if(verbose)
			cat("[parseValues] Parsing field", paste(.currentField), "\n")
		
		# convert values to the correct types
		.method <- .converters[[.currentField[[.sosParseFieldDefinition]]]]
		if(is.null(.method)) {
			# could still be a unit of measurement given, use as
			if(!is.na(.currentField[.sosParseFieldUOM])) {
				.method <- .converters[[.currentField[[.sosParseFieldUOM]]]]
				if(is.null(.method)) {
					# fallback option
					warning(paste("No converter for the unit of measurement ",
									.currentField[[.sosParseFieldUOM]],
									" with the definition ",
									.currentField[[.sosParseFieldDefinition]],
									"! Trying a default, but you can add one when creating a SOS using",
									"SosDataFieldConvertingFunctions()."))
					.method <- .converters[[.sosParseFieldUOM]]	
				}
			}
			else {
				warning(paste("No converter found for the given field",
								toString(.currentField)))
			}
		}
		
		if(is.null(.method)) {
			warning(paste("No converter found! Skipping field",
							as.character(fields[[.currentFieldIdx]]), "\n"))
			next;
		}
		
		if(verbose) {
			cat("[parseValues] Using converter function:\n")
			show(.method)
		}
				
		# do the conversion
		.currentValues <- .method(x = .currentValues, sos = sos)
		
		# bind new and existing data:
		if(verbose) cat("[parseValues] Binding additional data.frame for",
					.currentField[[.sosParseFieldName]],
					"-- value range", toString(range(.currentValues)), "\n")
		.newData <- data.frame(.currentValues)
		
		# create the names of the new data:
		.newDataName <- .currentField[[.sosParseFieldName]]		
		names(.newData) <- .cleanupColumnName(.newDataName)
		
		if(verbose) cat("[parseValues] Added column name:", names(.newData), "\n")

		# bind existing and new data column
		.data <- cbind(.data, .newData)
		
		if(verbose) {
			cat("[parseValues] The new bound data frame:\n")
			str(.data)
		}
		
		# add field information as attributes to the new column using human
		# readable names
		.addAttrs <- as.list(.currentField)
		names(.addAttrs) <- .sosParseFieldReadable[names(.currentField)]
		
		.lastColumn <- dim(.data)[[2]]
		.oldAttrs <- attributes(.data[,.lastColumn])
		
		attributes(.data[,.lastColumn]) <- c(as.list(.oldAttrs),
			.addAttrs)
		
		if(verbose) cat("[parseValues] Added attributes to new data:",
					toString(.addAttrs),
					"[ names: ", toString(names(.addAttrs)), "]",
					"\n[parseValues] Old attributes list is",
					toString(.oldAttrs),
					"\n[parseValues] New attributes list is",
					toString(attributes(.data[,.lastColumn])),
					"\n")
	}
	
	# remove id column
	if(verbose) cat("[parseValues] Removing temporary first column\n")
	.data <- .data[,!colnames(.data)%in%.tempId]
	
	if(verbose) cat("[parseValues] returning!\n")
	return(.data)
}

#
# Creates list of named character vectors with the information from swe:fields.
#
parseElementType <- function(obj, sos, verbose = FALSE) {
	.simpleDataRecord <- obj[[sweSimpleDataRecordName]]
	.dataRecord <- obj[[sweDataRecordName]]
	if(!is.null(.simpleDataRecord) || !is.null(.dataRecord)) {
		if(!is.null(.simpleDataRecord)) .dr <- .simpleDataRecord
		else .dr <- .dataRecord
		
		.fields <- .filterXmlChildren(node = .dr, childrenName = sweFieldName,
				includeNamed = TRUE)
		
		if(verbose) cat("[parseElementType] Got data record with",
					length(.fields), "fields. \n")
		
		# extract the fields, naming with attribute 'name'
		.parsedFields <- lapply(.fields, parseField, sos = sos,
				verbose = verbose)
		.names <- sapply(.parsedFields, "[", .sosParseFieldName)
		names(.parsedFields) <- .names
		
		if(verbose) cat("[parseElementType] Names of parsed fields:",
					names(.fields), "\n")
			
		return(.parsedFields)
	}
	else {
		stop(paste("Cannot parse swe:elementType, only children of type",
						sweSimpleDataRecordName, "and", sweDataRecordName,
						"are supported!"))
	}
}

#
# swe:encoding
#
parseEncoding <- function(obj, sos, verbose = FALSE) {
	.textBlock <- obj[[sweTextBlockName]]
	
	if(is.null(.textBlock)) {
		stop(paste("Cannot parse swe:encoding, only", sweTextBlockName,
						"is supported!"))
	}
	else {
		.tb <- parseTextBlock(.textBlock)
		return(.tb)
	}
}

################################################################################
# sub-parsing functions, not exchangeable via SosParsers

.sosParseFieldName <- ".sosParseFieldName"
.sosParseFieldDefinition <- ".sosParseFieldDefinition"
.sosParseFieldUOM <- ".sosParseFieldUOM"
.sosParseFieldCategoryName <- ".sosParseFieldCategoryName"
.sosParseFieldValue <- ".sosParseFieldValue"
.sosParseFieldCodeSpace <- ".sosParseFieldCodeSpace"

# human readable versions of the names:
.sosParseFieldReadable <- list(
		"name",
		"definition",
		"unit of measurement",
		"category name",
		"category value",
		"category code space")
names(.sosParseFieldReadable) <- list(
		.sosParseFieldName,
		.sosParseFieldDefinition,
		.sosParseFieldUOM,
		.sosParseFieldCategoryName,
		.sosParseFieldValue,#
		.sosParseFieldCodeSpace)

#
# Function creates a named character vector (using field names from variables
# like ".sosParseFieldXYZ") with all information stored in a swe:field.
#
parseField <- function(obj, sos, verbose = FALSE) {
	.field <- NULL
	.name <- xmlGetAttr(node = obj, name = "name")
	if(verbose) cat("[parseField] Parsing field description of ", .name, "\n")

	.noneText <- .filterXmlChildren(node = obj, childrenName = xmlTextNodeName,
			includeNamed = FALSE)
	.innerField <- .noneText[[1]]
	.innerFieldName <- xmlName(.innerField)
	
	# Available options: Time, Text, Quantity, Category
	# The parsed elements and fields are closely bound to 52N SOS (OMEncoder.java)
	if(.innerFieldName == sweTimeName) {
		.def <- xmlGetAttr(node = .innerField, name = "definition")
		
		.field <- c(.sosParseFieldName = .name, .sosParseFieldDefinition = .def)
	}
	else if (.innerFieldName == sweTextName) {
		.def <- xmlGetAttr(node = .innerField, name = "definition")
		
		.field <- c(.sosParseFieldName = .name, .sosParseFieldDefinition = .def)
	}
	else if (.innerFieldName == sweQuantityName) {
		.def <- xmlGetAttr(node = .innerField, name = "definition")
		
		if(!is.null(.innerField[[sweUomName]])) {
			.uom <- xmlGetAttr(node = .innerField[[sweUomName]], name = "code")
		}
		else {
			warning(paste("swe:Quantity given without unit of measurement:",
							.name))
			.uom <- NA_character_
		}
		.field <- c(.sosParseFieldName = .name, .sosParseFieldDefinition = .def,
				.sosParseFieldUOM = .uom)
	}
	else if (.innerFieldName == sweCategoryName) {
		.catName <- xmlGetAttr(node = .innerField, name = "name")
		.value <- xmlValue(.innerField[[sweValueName]])
		.codeSpace <- xmlGetAttr(node = .innerField[[sweCodeSpaceName]],
				name = "type")
		
		.field <- c(.sosParseFieldName = .name,
				.sosParseFieldCategoryName = .catName,
				.sosParseFieldValue = .value,
				.sosParseFieldCodeSpace = .codeSpace)
	}
	else if (.innerFieldName == sweBooleanName ||
			.innerFieldName == sweCountName) {
		warning("Parsing of the given swe:field", .innerFieldName,
				"is not implemented! Please extend the parsing funtion of swe:elementType.")
		.field <- c(name = .innerField)
	}
	else if (.innerFieldName == sweDataRecordName) {
		stop(paste("Parsing of nested swe:DataRecords is not supported!", 
						"Please extend the parsing funtion of swe:elementType."))
	}
	
	if(verbose) cat("[parseField] Parsed field", toString(.field), "\n")
	
	return(.field)
}

#
#
#
parseTextBlock <- function(obj) {
	.id <- xmlGetAttr(node = obj, name = "id", default = NA_character_)
	.tS <- xmlGetAttr(node = obj, name = "tokenSeparator")
	.bS <- xmlGetAttr(node = obj, name = "blockSeparator")
	.dS <- xmlGetAttr(node = obj, name = "decimalSeparator")
	
	.tb <- SweTextBlock(tokenSeparator = .tS, blockSeparator = .bS,
			decimalSeparator = .dS, id = .id)
	return(.tb)
}

#
#
#
parsePhenomenonProperty <- function(obj, sos, verbose = FALSE) {
	.obsProp <- NULL
	
	# check if reference or inline phenomenon
	.href <- xmlGetAttr(node = obj, name = "href")
	if(!is.null(.href)) {
		if(verbose) cat("[parsePhenomenonProperty] with reference", .href,
					"\n")
		
		.obsProp <- SwePhenomenonProperty(href = .href)
	}
	else {
		.noneText <- .filterXmlChildren(node = obj, xmlTextNodeName,
				includeNamed = FALSE)
		.compPhen <- .noneText[[1]]
		# 52N SOS only returns swe:CompositePhenomenon
		.name <- xmlName(.compPhen)
		if(verbose) cat("[parsePhenomenonProperty] inline with name", .name,
					"\n")
		
		if(.name == sweCompositePhenomenonName) {
			.phen <- parseCompositePhenomenon(.compPhen, sos = sos,
					verbose = verbose)
			.obsProp <- SwePhenomenonProperty(phenomenon = .phen)
		}
		else {
			warning(paste("[parsePhenomenonProperty] Unsupported observed property: ",
							.name, "\n"))
		}
	}
	
	return(.obsProp)
}

#
#
#
parseCompositePhenomenon <- function(obj, sos, verbose = FALSE) {
	.id <- xmlGetAttr(node = obj, name = "id", default = NA_character_)
	
	if(verbose) cat("[parseCompositePhenomenon] with id", .id, "\n")
	
	.dimension <- as.integer(
			xmlGetAttr(node = obj, name = "dimension", default = NA_character_))
	.name <- xmlValue(obj[[gmlNameName]])
	
	.components <- lapply(obj[sweComponentName], parseComponent, sos = sos,
			verbose = verbose)
	
	if(verbose) cat("[parseCompositePhenomenon]", length(.components),
				"components parsed.\n")
	
	# optional:
	.description <- NA_character_
	if(!is.null(obj[[gmlDescriptionName]])) {
		.description <- parsePhenomenonProperty(obj[[sweBaseName]])
	}
	.base <- NULL
	if(!is.null(obj[[sweBaseName]])) {
		.base <- parsePhenomenonProperty(obj[[sweBaseName]])
	}
	
	.compPhen <- SweCompositePhenomenon(id = .id, name = .name, 
			description = .description, dimension = .dimension,
			components = .components, base = .base)
	
	return(.compPhen)
}

#
#
#
parseComponent <- function(obj, sos, verbose = FALSE) {
	# 52N SOS only sets the href property on swe components, but still reuse function
	.component <- parsePhenomenonProperty(obj)
	return(.component)
}

#
#
#
parseSwePosition <- function(obj, sos, verbose = FALSE) {
	.rF <- xmlGetAttr(node = obj, name = "referenceFrame")
	if(verbose) cat("[parseSwePosition] with referenceFrame", .rF, "\n")
	
	.location <- obj[[sweLocationName]]
	.parser <- sosParsers(sos)[[sweLocationName]]
	
	.pos <- .parser(.location, sos = sos, verbose = verbose)
	
	.oldAttrs <- attributes(.pos)
	attributes(.pos) <- c(.oldAttrs, list(referenceFrame = .rF))
	
	return(.pos)
}

#
#
#
parseLocation <- function(obj, sos, verbose = FALSE) {
	.vector <- obj[[sweVectorName]]
	.id <- xmlGetAttr(node = obj, name = "id")
	if(verbose) cat("[parseLocation] with id", .id, "\n")
	
	.parser <- sosParsers(sos)[[sweVectorName]]
	
	.pos <- .parser(.vector, sos = sos, verbose = verbose)
	return(.pos)
}

#
#
#
parseVector <- function(obj, sos, verbose = FALSE) {
	.children <- .filterXmlChildren(node = obj, childrenName = sweCoordinateName)
	
	.coords <- list()
	.names <- list()
	for (i in seq(1, length(.children))) {
		.c <- .children[[i]]
		
		.coord <- parseCoordinate(.c, sos = sos, verbose = verbose)
		.coords <- c(.coords, list(.coord))
		.names <- c(.names, list(.coord[["axisID"]]))
		
		if(verbose) cat("[parseVector] parsed coordinate: ", toString(.coord),
					"\n")
	}
	
	names(.coords) <- .names

	return(.coords)
}

#
#
#
parseCoordinate <- function(obj, sos, verbose = FALSE) {
	.name <- xmlGetAttr(node = obj, name = "name", default = NA_character_)
	if(verbose) cat("[parseCoordinate] with name", .name, "\n")

	.quantity <- obj[[sweQuantityName]]
	.axisID <- xmlGetAttr(.quantity, name = "axisID", default = NA_character_)
	if(verbose) cat("[parseCoordinate] axisID: ", .axisID, "\n")
	
	.uomCode <- NA_character_
	if(!is.null(.quantity[[sweUomName]]))
		.uomCode <- xmlGetAttr(node = .quantity[[sweUomName]], name = "code",
				default = NA_character_)
	if(verbose) cat("[parseCoordinate] uomCode: ", .uomCode, "\n")
	
	.value <- NA_character_
	if(!is.null(.quantity[[sweValueName]]))
		.value <- as.double(xmlValue(.quantity[[sweValueName]]))
	if(verbose) cat("[parseCoordinate] value: ", .value, "\n")
	
	return(list(name = .name, axisID = .axisID, uomCode = .uomCode,
					value = .value))
}
