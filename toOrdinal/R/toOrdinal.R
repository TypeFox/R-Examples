`toOrdinal` <-
function(
	cardinal_number,
	language="English",
	convert_to="ordinal_number") {


	### Utility function

	strtail <- function(s, n=1) {
		if(n < 0) substring(s, 1-n)
		else substring(s, nchar(s)-n+1)
	}


	### Argument tests

	supported_languages_ordinal_number <- c("ENGLISH", "FRENCH", "GERMAN", "SPANISH", "SWEDISH")
	supported_languages_ordinal_word <- ""
	if (floor(cardinal_number)!=cardinal_number | cardinal_number < 0) stop("Number supplied to 'toOrdinal' must be a positive integer.", call.=FALSE)


	#######################################################
	###
	### convert_to ordinal_number
	###
	#######################################################

	if (identical(toupper(convert_to), "ORDINAL_NUMBER")) {

		if (!toupper(language) %in% supported_languages_ordinal_number) stop(paste("Language supplied (", language, ") is currently not supported by toOrdinal for conversion to an 'ordinal_number'. Currently supported languages include: ", paste(supported_languages_ordinal_number, collapse=", "), ". Please submit pull requests to https://github.com/CenterForAssessment/toOrdinal/pulls for additional language support.", sep=""), call.=FALSE)


		### ENGLISH

		if (toupper(language)=="ENGLISH") {
			tmp <- strtail(as.character(cardinal_number), 2)
			if (tmp %in% c('1', paste(c(0, 2:9), 1, sep=""))) tmp.suffix <- "st"
			if (tmp %in% c('2', paste(c(0, 2:9), 2, sep=""))) tmp.suffix <- "nd"
			if (tmp %in% c('3', paste(c(0, 2:9), 3, sep=""))) tmp.suffix <- "rd"
			if (tmp %in% c('11', '12', '13')) tmp.suffix <- "th"
			if (tmp %in% c('4', paste(0:9, 4, sep=""))) tmp.suffix <- "th"
			if (tmp %in% c('5', paste(0:9, 5, sep=""))) tmp.suffix <- "th"
			if (tmp %in% c('6', paste(0:9, 6, sep=""))) tmp.suffix <- "th"
			if (tmp %in% c('7', paste(0:9, 7, sep=""))) tmp.suffix <- "th"
			if (tmp %in% c('8', paste(0:9, 8, sep=""))) tmp.suffix <- "th"
			if (tmp %in% c('9', paste(0:9, 9, sep=""))) tmp.suffix <- "th"
			if (tmp %in% c('0', paste(0:9, 0, sep=""))) tmp.suffix <- "th"
		}


		### FRENCH

		if (toupper(language)=="FRENCH") {
			if (cardinal_number==1) tmp.suffix <- "re" else tmp.suffix <- "e"
		}


		### GERMAN

		if (toupper(language)=="GERMAN") {
			if (cardinal_number >=1 & cardinal_number <= 19) tmp.suffix <- "te"
			if (cardinal_number >= 20) tmp.suffix <- "ste"
			}


		### SPANISH

		if (toupper(language)=="SPANISH") {
			tmp <- strtail(as.character(cardinal_number), 1)
			if (tmp %in% c('1', '3')) tmp.suffix <- ".er"
			if (tmp %in% c('0', '2', '4', '5', '6', '7', '8', '9')) tmp.suffix <- ".\u00BA"
		}


		### SWEDISH

		if (toupper(language)=="SWEDISH") {
		 	tmp_1char <- strtail(as.character(cardinal_number), 1)
		 	tmp_2char <- strtail(as.character(cardinal_number), 2)
		 	if (tmp_1char %in% c('0', '3', '4', '5', '6', '7', '8', '9') | tmp_2char %in% c('11', '12')) {
				tmp.suffix <- ":e"
			} else if (tmp_1char %in% c('1', '2')) {
				tmp.suffix <- ":a"
			}
		}


		### TURKISH

		if (toupper(language)=="TURKISH") {
		}

		return(paste(cardinal_number, tmp.suffix, sep=""))

	} ### if (identical(toupper(convert_to), "ORDINAL_NUMBER"))


	######################################################################
	###
	### convert_to ordinal_word
	###
	######################################################################

	if (identical(toupper(convert_to), "ORDINAL_WORD")) {

		if (!toupper(language) %in% supported_languages_ordinal_word) stop(paste("Language supplied (", language, ") is currently not supported by toOrdinal for conversion to an 'ordinal_word'. Currently supported languages include: ", paste(supported_languages_ordinal_word, collapse=", "), ". Please submit pull requests to https://github.com/CenterForAssessment/toOrdinal/pulls for additional language support.", sep=""), call.=FALSE)


		### ENGLISH




	} ### if (identical(toupper(convert_to), "ORDINAL_WORD"))
} ### END toOrdinal
