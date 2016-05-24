# SETMETHOD CONVERT.UNITS####################################################3
############# convert.units: a function that converts units from one to another, if possible.
#    x         = a numeric vector with values to be converted
#    tounit    = a character vector of the new units to be used. Must be found from the To column from the table in [[Unit conversions]].
#    fromunit  = a character vector or factor with the current units
convert.units <- function(x, tounit = c("kg", "s", "m", "m3", "J", "W", "A", "V", "C", "N", "Pa", "Hz", "mol"), fromunit = NULL) {
	if(is.null(tounit)) {tounit <- c("kg", "s", "m", "m3", "J", "W", "A", "V", "C", "N", "Pa", "Hz", "mol")}
	if(is.null(fromunit)) { # Do nothing if from-units are not defined.
		out <- data.frame(Unit = fromunit, Result = x)
	} else {
		conversions <- tidy(opbase.data("Op_en5475")) # Set up a full unit conversion table.
		conversions$Result <- suppressWarnings(as.numeric(as.character(conversions$Result)))
		#colnames(conversions)[colnames(conversions) == "type"] <- "Type"
		prefixes <- conversions[conversions$Type == "Prefix", ] # Combine all prefixes with all units in From column.
		colnames(prefixes) <- paste("Prefix.", colnames(prefixes), sep = "")
		conversions <- merge(prefixes, conversions[conversions$Type == "Unit", ])
		conversions$From <- paste(conversions$Prefix.From, conversions$From, sep = "")
		conversions$Result <- conversions$Result * conversions$Prefix.Result
		conversions <- conversions[,c("From", "To", "Result")]
		conversions <- merge(conversions, conversions, by = "To") # Create all possible from-to pairs.
		conversions$Result <- conversions$Result.x / conversions$Result.y
		conversions <- unique(conversions[c("From.x", "From.y", "Result")])
		colnames(conversions) <- c("From", "To", "Result")
		coefficients <- data.frame()
		for(i in levels(as.factor(fromunit))) { # Look through each different "From" unit in the data.
			coefficient <- 1
			outto <- ""
			for(j in strsplit(i, split = " ")[[1]]) { # Look though each part of a composite unit.
				if(gsub("^/", "", j) != j) { # If unit is in denominator, use inverse.
					j <- gsub("^/", "", j)
					exponent <- -1
					slash <- "/"
				} else {
					exponent <- 1
					slash <- ""
				}
				conversions.j <- merge(conversions, j, by.x = "From", by.y = "y")
				conversions.j <- conversions.j[conversions.j$To %in% tounit, ]
				if(nrow(conversions.j) > 0) {
					coefficient <- coefficient * conversions.j$Result^exponent
					outto <- paste(outto, " ", slash, conversions.j$To, sep = "")
				}
			}
			coefficients <- rbind(coefficients, data.frame(
				From = i,
				To = gsub("^ ", "", outto), 
				Result = coefficient
			))
		}
		# Combine coefficients with data.
		out <- merge(fromunit, coefficients, by.x = "x", by.y = "From", sort = FALSE)
		out <- data.frame(Unit = out$To, Result = x * out$Result)
	}
	return(out)
}

setGeneric("convert.units")

setMethod(
	f = "convert.units",
	signature = signature(x = "ovariable"),
	definition = function(x, tounit)
	{
		if("Unit" %in% colnames(x@output) & !is.null(tounit)) {
			x@output[c("Unit", "Result")] <- callGeneric(x = x@output$Result, tounit = tounit, fromunit = x@output$Unit)#, expo = expo)
		}
		return(x)
	}
)

setMethod(
	f = "convert.units",
	signature = signature(x = "data.frame"),
	definition = function(x, tounit)
	{
		if("Unit" %in% colnames(x) & "Result" %in% colnames(x) & !is.null(tounit)) {
			x[c("Unit", "Result")] <- callGeneric(x = x$Result, tounit = tounit, fromunit = x$Unit)#, expo = expo)
		}
		return(x)
	}
)