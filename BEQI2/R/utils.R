#' 	Remove Redundant Spaces
#'
#' 	This function removes redundant spaces from character vectors
#'
#' 	@param x character vector
#'
#'  @return  character vector without trailing or multiple spaces
#'  @examples stopifnot(BEQI2:::stripSpaces(" Hello  World  ") == "Hello World")
stripSpaces <- 
function(x) {
	x <- gsub(pattern = " {2,}",   replacement = " ", x = x)
	x <- gsub(pattern = "^ +| +$", replacement = "",  x = x)
	x
}


#'  Test for Azoic Samples
#'
#' 	Case-insensitive test for taxa starting with 'azoi'
#'
#' 	@param x character vector containing taxa
#'
#'  @return  logical vector, with elements \code{TRUE} for azoic samples, 
#'      and \code{FALSE} otherwise.
isAzoic <- 
function(x) {
    grepl(pattern = "^Azoi", x = x, ignore.case = TRUE) 
}




#' 	Read BEQI Settings File
#'
#' 	This function reads BEQI settings files (JSON)
#'
#' 	@param filename name of BEQI input file (\code{character})
#'
#' 	@details The function performs the following tasks:
#' 	\itemize{
#' 		\item{checks the existence of \code{filename};}
#'  	\item{reads JSON file while ignoring C-style comments;}
#'      \item{checks avaiability of required keys in the JSON-file}
#'  	\item{checks values in JSON-file}
#'  }
#'
#'	@import jsonlite
#'    
#'  @export
readSettings <- 
function(filename) {

    # check existence of settings file
	if (!file.exists(filename)) {
		stop("File not found", call. = FALSE)
	}

	# read settings file
	settings <- readLines(con = filename, warn = FALSE)

	# remove all C-style comments (//)
	# Note: comments are formally not part of the JSON specification.
	settings <- sub(pattern = "//.*$", replacement = "", x = settings)

	# parse JSON
	if (!validate(settings)) {
		stop(
			sprintf(
				"Errors found in %s. Check JSON-format (e.g. brackets, braces, trailing comma's)", 
				sQuote(filename)
			),  
			call. = FALSE
		)
	}
	settings <- fromJSON(settings)
    
    # check if required keys are available
    requiredKeys <- c("title", "user", "date", "files",
                      "months", "pooling", "genusToSpeciesConversion")
    names(settings) <- tolower(names(settings))
    found <- tolower(requiredKeys) %in% names(settings)
    if (any(!found)) {
        stop(
            sprintf(
                fmt = "key %s is missing in JSON-file %s\n(see package vignette)", 
                toString(sQuote(requiredKeys[!found])), 
                sQuote(basename(filename))
            ), 
            call. = FALSE
        )
    }
    
    # check required files
    requiredKeys <- c("BEQI2", "SpeciesNames", "Ecotopes")
    names(settings$files) <- tolower(names(settings$files))
    found <- tolower(requiredKeys) %in% names(settings$files)
    if (any(!found)) {
        stop(
            sprintf(
                fmt = "key files:%s is missing in JSON-file %s\n(see package vignette)", 
                toString(sQuote(requiredKeys[!found])), 
                sQuote(basename(filename))
            ), 
            call. = FALSE
        )
    }

    # check months
    if (!is.integer(settings$months) | length(settings$months) != 2L) {
    	stop(
            "key 'months' should be an integer vector of length 2", 
             call. = TRUE
        )
    }
	if (!all(settings$months %in% 1:12)) {
		stop("elements of key 'months' should be in [1, 12]", call. = TRUE)
	}
    if ((settings$months[2] - settings$months[1]) < 0) {
		stop(
            "First month to analyse should be smaller than or equal to last month", 
            call. = TRUE
        )
	}

    # check data pooling
    names(settings$pooling) <- tolower(names(settings$pooling))
    if (!is.logical(settings$pooling$enabled)) {
        stop(
            "key 'pooling:enabled' should be either 'true' or 'false'", 
            call. = TRUE
        )
    }
    if (settings$pooling$enabled) {
        settings$pooling$targetarea <- range(settings$pooling$targetarea)
        if (!is.numeric(settings$pooling$targetarea) | 
            (length(settings$pooling$targetarea) != 2L)) {
        	stop(
                "key 'pooling:targetArea' should be a numeric vector of length 2", 
                call. = TRUE
            )
        }
        settings$pooling$randomseed <- as.integer(settings$pooling$randomseed)
        if (!is.integer(settings$pooling$randomseed)) {
            stop(
                "key 'pooling:randomSeed' needs to be an integer vector of length 1", 
                call. = TRUE
            )
        }
    }

    # genus to species conversion
    if (!is.logical(settings$genustospeciesconversion)) {
        stop(
            "key 'genusToSpeciesConversion' should be either 'true' or 'false'", 
            call. = TRUE
        )
    }
    
	# return results
	settings
}


#' 	Read BEQI input files
#'
#' 	This function reads and checks BEQI input files. The format has been
#'	specified in Van Loon (2013).
#'
#' 	@param filename name of BEQI input file (character)
#'
#' @details The function performs the following tasks:
#' 	\itemize{
#' 		\item{checks the existence of \code{filename};}
#'  	\item{checks availablitity of required columns (case insensitive);}
#'      \item{make column names with aggregation data case-insensitive;}
#'  	\item{removes redundant spaces;}
#'      \item{checks if DATE-field adheres to ISO 8601 (YYYY-mm-dd);}
#'  	\item{constructs a unique identifier \code{ID} by concatenating 
#'          columns \code{OBJECTID} and \code{DATE};}
#'      \item{aggregate (by summation) VALUE-fields of records that only differ 
#'          in VALUE-field value;}
#'      \item{checks that each \code{ID} has a unique \code{AREA};}
#'      \item{checks azoic samples for VALUE=0;}
#'      \item{removes records with VALUE=0, not belonging to azoic samples;}
#'      \item{checks VALUE-field on missing values;}
#'      \item{checks if VALUE-field is an integer;}
#'  }
#'
#' 	@references Willem van Loon, 2013. BEQI2 INPUT FORMAT
#'
#' 	@export
readBEQI <-
function(filename) {

	# check if 'filename' exists
	if (!file.exists(filename)) {
		stop(
			sprintf(fmt = "File %s not found", sQuote(filename)),
			call. = FALSE
		)
	}

	# read file 'filename'
	d <- try(read.csv(file = filename, as.is = TRUE), silent = TRUE)
	if (inherits(d, "try-error")) {
		stop(
			sprintf(fmt = "Errors occurred while reading %s", sQuote(filename)),
			call. = FALSE
		)
	}

	# check column names (case insensitive)
	requiredColumns <- c("OBJECTID", "ECOTOPE", "SAMPLEID", "TAXON", 
                         "CHAR", "SAMPDEV", "AREA", "DATE", "VALUE")
	missingColumns <- setdiff(requiredColumns, toupper(names(d)))
	if (length(missingColumns) > 0L) {
		stop(
            sprintf(
                fmt = "The following columns are missing in file\n%s: %s", 
                sQuote(basename(filename)),
                toString(missingColumns)
            ),
            call. = FALSE
        )
	}
	names(d) <- toupper(names(d))
    
    # check VALUE field on missing values
    index <- which(is.na(d$VALUE))
    if (length(index) != 0L) {
        stop(
            sprintf(
                fmt = paste0(
                    "The value in field %s is missing for %i records.\n",
                    "The record indices are:\n%s"
                ),
                sQuote("VALUE"),
                length(index),
                toString(index)
            ), 
            call. = FALSE
        )
    }
    
    # check VALUE field on integers
    d$VALUE <- suppressWarnings(as.numeric(d$VALUE))
    isInteger <- (abs(round(d$VALUE) - d$VALUE) < .Machine$double.eps)
    index <- which(is.na(isInteger) | !isInteger)
    if (length(index) != 0L) {
        stop(
            sprintf(
                fmt = paste0(
                    "Field %s contains %i non-integer values.\n",
                    "The record indices are:\n%s"
                ),
                sQuote("VALUE"),
                length(index),
                toString(index)
            ), 
            call. = FALSE
        )
    }
    
	# remove redundant spaces
	d <- as.data.frame(
		lapply(X = d, FUN = function(x) {
			if (is.character(x)) {
				x <- stripSpaces(x)
			}
			x
		}),
        stringsAsFactors = FALSE
	)
    
    # harmonize columns to make case-insensitive matching possible
    d <- as.data.frame(
		lapply(X = d, FUN = function(x) {
			if (is.character(x)) {
				x <- harmonize(x)
			}
			x
		}),
        stringsAsFactors = FALSE
	)
    
	# coerce date to class 'Date'
	res <- try(as.Date(d$DATE), silent = TRUE)
	if (inherits(res, "try-error") | any(is.na(res))) {
		stop(
			"Invalid date formats found. Please adhere to ISO 8601 (YYYY-mm-dd)",
			call. = FALSE
		)
	}
	d$DATE <- res

    # add (unique) identifier field
    d$ID <- paste(d$OBJECTID, d$SAMPLEID, d$DATE, sep = "/")

    # check if areas are unique for a specific ID
    n <- tapply(
        X = d$AREA, 
        INDEX = d$ID, 
        FUN = function(x){length(unique(x))}
    )
    if (!isTRUE(all(n == 1L))) {
        n <- names(n)[n > 1L]
        warning(
            sprintf(
                fmt = paste0(
                    "AREA-field is not unique for %i samples. ",
                    "These will be removed.\n",
                    "These samples have the following identifiers:\n%s"
                ),
                length(n),
                toString(sQuote(n))
            ), 
            call. = FALSE
        )
        d <- d[!(d$ID %in% n), ]
    }

    # aggregate (by summation) the VALUE-fields of records that at most only 
    # differ in VALUE-field value. These samples are the result of duplo/triplo 
    # samples and need to be combined.
    # (NB: the sample areas should not be summed!)
    VALUE <- d$VALUE
    d <- subset(x = d, select = -VALUE)
    d$INDEX <- apply(X = d, MARGIN = 1L, FUN = paste0, collapse = "")
    VALUE <- tapply(X = VALUE, INDEX = d$INDEX, FUN = sum)
    d <- unique(d)
    d$VALUE <- as.integer(VALUE)[match(x = d$INDEX, table = names(VALUE))]
    d$INDEX <- NULL

    # handle azoic samples
    index <- which(isAzoic(x = d$TAXON))
    if (length(index) != 0L) {
        if (isTRUE(any(d$VALUE[index] != 0L))) {
            warning(
                paste(
                    "Azoic sample(s) found with non-zero VALUE-field.",
                    "These abundances will be set to zero"
                ),
                call. = FALSE
            )
        }
        d$VALUE[index] <- 0L
    }
    
    # remove non-Azoic records with zero counts as these are redundant 
    # and won't affect the results
    isZero <- (d$VALUE == 0L) & !isAzoic(x = d$TAXON)
    if (any(isZero)) {
        warning(
            sprintf(
                paste(
                    "Non-azoic records (n=%s) found with zero VALUE-field.",
                    "These records are redundant and will be excluded."
                ),
                sum(isZero)
            ),
            call. = FALSE
        )
        d <- d[!isZero, ]
    }
    
    # final check
    if (nrow(d) == 0L) {
        stop(
            sprintf(
                fmt = "No valid records found in: %s",
                sQuote(filename)
            ), 
            call. = FALSE
        )
    }

	# return content of file
	d
}







#' 	Read TWN Data
#'
#' 	This function reads files in the Taxa Waterbeheer Nederland (TWN) format.
#'
#' 	@details The function adds a new column \code{taxon}. Its contents depending on TWN-status:
#' 	\itemize{
#' 		\item{status = 10} {taxonname}
#'  	\item{status = 20} {prefername}
#'  	\item{status = 80} {parentname}
#' 	}
#'
#' 	@param filename name of TWN  listBEQI input file (\code{character})
#'
#'  @return a \code{data.frame} with four columns: taxonname, taxongroup, 
#'      taxonlevel, taxon
#'
#' 	@references \url{sofus.ecosys.nl/taxabase.htm}
#'  @references \url{www.aquo.nl/faq/faq-twn}
#'
#' 	@export
readTWN <-
function(filename) {

	# check if 'filename' exists
	if (!file.exists(filename)) {
		stop(
			sprintf(fmt = "File %s not found", sQuote(filename)),
			call. = FALSE
		)
	}

	# read file 'filename'
	d <- try(read.csv(file = filename, as.is = TRUE), silent = TRUE)
	if (inherits(d, "try-error")) {
		stop(
			sprintf(fmt = "Errors occurred while reading %s", sQuote(filename)),
			call. = FALSE
		)
	}

	# check column names
	requiredColumns <- c("status", "taxonname", "taxongroup", "prefername", 
                         "parentname", "taxonlevel")
	missingColumns <- setdiff(requiredColumns, tolower(names(d)))
	if (length(missingColumns) > 0L) {
		stop(
            sprintf(
                fmt = "The following columns are missing: %s", 
                toString(missingColumns)
            ),
            call. = FALSE
        )
	}
	names(d) <- tolower(names(d))

	# select only columns of interest
	d <- d[, requiredColumns]

	# keep only status codes:
	#	10: preferred name 
	#	20: synonym
	#	80: non-taxonomic species group
	# see also www.aquo.nl/faq/faq-twn
	d <- d[d$status %in% c(10L, 20L, 80L), ]

	# remove redundant spaces 
	# (including leading and trailing spaces)
	d$taxonname  <- stripSpaces(d$taxonname)
	d$prefername <- stripSpaces(d$prefername)
	d$parentname <- stripSpaces(d$parentname)
    d$taxongroup <- stripSpaces(d$taxongroup)

	# construct taxon
	d$taxon <- NA_character_
	d$taxon[d$status == 10L] <- d$taxonname [d$status == 10L]
	d$taxon[d$status == 20L] <- d$prefername[d$status == 20L]
	d$taxon[d$status == 80L] <- d$parentname[d$status == 80L]

	if (any(is.na(d$taxon))) {
		warning(
			sprintf(
				"A total of %s taxon names cannot be converted", 
				sum(is.na(d$taxon))
			), 
			call. = TRUE
		)
	}

	# create ordered factor of taxon levels
	d$taxonlevel <- factor(
		x = d$taxonlevel, 
		levels = c(
			"Regio", 
			"Regnum", 
			"Phylum", "Subphylum", 
			"Classis", "Subclassis", "Infraclassis", 
			"Ordo", "Subordo", "Infraordo", 
			"Superfamilia",	"Familia", "Subfamilia", 
			"Tribe", 
			"Genus", "Genus combi", "Subgenus", 
			"Species", "Species combi", "Subspecies", 
			"Varietas", 
			"Forma"),
		ordered = TRUE
	)

	# return subset
	d[, c("taxonname", "taxongroup", "taxonlevel", "taxon")]
}






#'  Read AMBI Sensitivity Data
#'
#' 	This function reads and checks files with AMBI sensitivity data. If 
#'     \code{filename = NULL} Borja's data will be read
#'
#' 	@param filename name of the AMBI sensitivity file (character)
#'
#' @details The function performs the following tasks:
#' 	\itemize{
#' 		\item{checks the existence of \code{filename};}
#'  	\item{checks availablitity of required columns (case insensitive);}
#'  	\item{removes redundant spaces;}
#'  	\item{removes duplicated records.}
#'  }
#'
#' 	@export
readAMBI <-
function(filename = NULL) {

    # check 'filename'
    if (is.null(filename)) {
        d <- readRDS(file = system.file("extdata/AMBI.rds", package = "BEQI2"))
    } else {

        if (!(file.exists(filename))) {
    		stop(
    			sprintf(fmt = "File %s not found", sQuote(filename)),
    			call. = FALSE
    		)
    	}

    	d <- try(read.csv(file = filename, as.is = TRUE), silent = TRUE)
    	if (inherits(d, "try-error")) {
    		stop(
    			sprintf(fmt = "Errors occurred while reading %s", sQuote(filename)),
    			call. = FALSE
    		)
    	}
    }

    # check column names (case insensitive)
    names(d) <- toupper(names(d))
	requiredColumns <- "SUBMITTED.NAME"
	missingColumns <- setdiff(requiredColumns, names(d))
	if (length(missingColumns) > 0L) {
		stop(
            sprintf(
                fmt = "The following column is missing: %s", 
                toString(missingColumns)
            ),
            call. = FALSE
        )
	}
    
    # check number of columns 
    # (should be 2 since the name of the second column is not fixed)
    if (ncol(d) != 2L) {
        stop(
            "Two columns expected 'SUBMITTED.NAME' and 'AMBI'",
            call. = FALSE
        )
    }

	# remove redundant spaces
	d <- data.frame(
		lapply(
            X = d, 
            FUN = function(x) {
    			if (is.character(x)) {
    				x <- stripSpaces(x)
    			}
    			x
    		}
        ), 
        stringsAsFactors = FALSE
	)

	# remove duplicated records
	if (anyDuplicated(d)) {
		isDuplicated <- duplicated(d)
		message(sprintf(fmt = "Number of duplicated records: %i", sum(isDuplicated)))
		message("These will be removed")
		d <- d[!isDuplicated, ]
	}
    
    # check on duplicated taxa
    if (anyDuplicated(d$SUBMITTED.NAME)) {
    	index <- which(duplicated(d$SUBMITTED.NAME))
		stop(
            sprintf(
                fmt = "Duplicated taxonnames found: %s", 
                toString(d$SUBMITTED.NAME[index])
            ),
            call. = FALSE
		)
    }

    # return content of file
	d
}



#'  Read Infaunal Trophic Index File
#'
#'  This function reads and checks files containing Infaunal trophic index 
#'  data (Gittenberger et al., 2011)
#'
#' 	@param filename name of the ITI file (character)
#'
#'  @details The function performs the following tasks:
#' 	\itemize{
#' 		\item{checks the existence of \code{filename};}
#'  	\item{checks availablitity of required columns (case insensitive);}
#'  	\item{removes redundant spaces;}
#'  	\item{removes duplicated records.}
#'      \item{checks if all ITI classes are I, II, III, or IV}
#'  }
#'
#' 	@export
readITI <-
function(filename) {
    read_ITI_FIBI(filename = filename, which = "ITI")
}


#'  Read Freshwater Inflow Biotic Index File
#'
#'
#' @param filename name of the FIBI file (character)
#'
#' @details The function performs the following tasks:
#' 	\itemize{
#' 		\item{checks the existence of \code{filename};}
#'  	\item{checks availablitity of required columns (case insensitive);}
#'  	\item{removes redundant spaces;}
#'  	\item{removes duplicated records.}
#'      \item{checks if all FIBI classes are I, II, III, or IV}
#'  }
#'  
#' 	@export
readFIBI <-
function(filename) {
    read_ITI_FIBI(filename = filename, which = "FIBI")
}


read_ITI_FIBI <- 
function(filename, which = c("ITI", "FIBI")) {
    
    # check argument 'which'
    which <- match.arg(which)

	# check if 'filename' exists
	if (!file.exists(filename)) {
		stop(
			sprintf(fmt = "File %s not found", sQuote(filename)),
			call. = FALSE
		)
	}

	# read file 'filename'
	d <- try(read.csv(file = filename, as.is = TRUE), silent = TRUE)
	if (inherits(d, "try-error")) {
		stop(
			sprintf(fmt = "Errors occurred while reading %s", sQuote(filename)),
			call. = FALSE
		)
	}

	# check column names (case insensitive)
    names(d) <- toupper(names(d))
	requiredColumns <- c("SUBMITTED.NAME", which)
	missingColumns <- setdiff(requiredColumns, names(d))
	if (length(missingColumns) > 0L) {
		stop(
            sprintf(
                fmt = "The following columns are missing: %s", 
                toString(missingColumns)
            ),
            call. = FALSE
        )
	}

	# remove redundant spaces
	d <- data.frame(
		lapply(
            X = d, 
            FUN = function(x) {
    			if (is.character(x)) {
    				x <- stripSpaces(x)
    			}
    			x
    		}
        ), 
        stringsAsFactors = FALSE
	)

	# remove duplicated records
	if (anyDuplicated(d)) {
		isDuplicated <- duplicated(d)
		message(sprintf(fmt = "Number of duplicated records: %i", sum(isDuplicated)))
		message("These will be removed")
		d <- d[!isDuplicated, ]
	}

    # check on duplicated taxa
    if (anyDuplicated(d$SUBMITTED.NAME)) {
    	index <- which(duplicated(d$SUBMITTED.NAME))
		stop(
            sprintf(
                fmt = "Duplicated taxonnames found: %s", 
                toString(d$SUBMITTED.NAME[index])
            ),
            call. = FALSE
		)
    }

    # check if all classes are in {I, II, III, IV}
    expectedClasses <- c("I", "II", "III", "IV")
    if (!all(d[[which]] %in% expectedClasses)) {
		stop(
            sprintf(
                fmt = "%s classes should be in {%s}", 
                sQuote(which),
                toString(expectedClasses)
            ),
            call. = FALSE
		)
    }

    # return content of file
	d
}


#'  Read Ecotopes References File
#'
#'  This function reads and checks files with reference values
#'
#' 	    @param filename name of the ecotopes reference file (character)
#'      @param extra additional user-defined indices to be checked (\code{character}, see details)
#'
#' @details The function performs the following tasks:
#' 	\itemize{
#' 		\item{checks the existence of \code{filename};}
#'  	\item{checks availablitity of required columns (case insensitive);}
#'  	\item{removes redundant spaces}
#'  	\item{removes duplicated records}
#'  }
#'  
#'  Argument \code{extra} is a \code{character} vector of additional benthic 
#'  indices to be checked for. For example, if \code{extra = "ITI"}, then
#'  the ecotope reference file should also contain columns ITIREF and ITIBAD.
#'  
#'  The format of the ecotopes reference file is documented in the 
#'  BEQI2-package vignette.
#'
#'  @references Van Loon, W, 2013. Loon2013-BEQI2-Specs-Ecotopes-27nov.doc
#'  
#' 	@export
readERF <-
function(filename, extra = NULL) {

	# check if 'filename' exists
	if (!file.exists(filename)) {
		stop(
			sprintf(fmt = "File %s not found", sQuote(filename)),
			call. = FALSE
		)
	}

	# read file 'filename'
	d <- try(
        read.csv(file = filename, as.is = TRUE), 
        silent = TRUE
    )
	if (inherits(d, "try-error")) {
		stop(
			sprintf(fmt = "Errors occurred while reading %s", sQuote(filename)),
			call. = FALSE
		)
	}

	# check column names (case insensitive)
	requiredColumns <- c("OBJECTID", "ECOTOPE", "RELAREA", 
        "SREF", "SBAD", "HREF", "HBAD", "AMBIREF", "AMBIBAD")
    if (!is.null(extra)) {
        requiredColumns <- c(requiredColumns,
            paste0(rep(extra, each = 2), c("REF", "BAD"))
        )
    }
	missingColumns <- setdiff(requiredColumns, names(d))
	if (length(missingColumns) > 0L) {
		stop(
            sprintf(
                fmt = "The following columns are missing: %s", 
                toString(missingColumns)
            ),
            call. = FALSE
        )
	}

    # remove redundant spaces
	d <- data.frame(
		lapply(
            X = d, 
            FUN = function(x) {
    			if (is.character(x)) {
    				x <- stripSpaces(x)
    			}
    			x
    		}
        ), 
        check.names = FALSE,
        stringsAsFactors = FALSE
	)

	# remove duplicated records
	if (anyDuplicated(d)) {
		isDuplicated <- duplicated(d)
		message(sprintf(fmt = "Number of duplicated records: %i", sum(isDuplicated)))
		message("These will be removed")
		d <- d[!isDuplicated, ]
	}

    # harmonize columns to make case-insensitive matching possible
    d <- as.data.frame(
    	lapply(X = d, FUN = function(x) {
			if (is.character(x)) {
				x <- harmonize(x)
			}
			x
		}),
        stringsAsFactors = FALSE
	)

    # return content of file
	d
}









#'  Test if a Value is in an Interval
#'
#' 	This function tests if values are part of a closed interval.
#'
#' 	@param e1 numeric value
#'  @param e2 numeric interval
#'
#'  @return  TRUE if the interval includes the value, FALSE otherwise
#'  @examples \dontrun{3 %inInterval% c(1, 4)}
"%inInterval%" <- function(e1, e2) {
    e2 <- range(e2)
    e1 >= e2[1] & e1 <= e2[2]
}


#'  Pooling
#'
#'  This function randomly assigns samples to pools of approximately equal area
#'
#' 	@param sampleId identifier
#'  @param area sample corresponding to \code{sampleId}
#'  @param targetArea vector of length 2 containing the lower and upper bound 
#'      of the pooled area
#'  @param maxTry maximum number of unsuccessful pooling tries before the 
#'      algorithm gives up.
#'
#'  @return vector with idenitifiers (integers) indicating the pool to which 
#'      each sample belongs (NA for samples that could not be pooled)
#'      
#'  @export
pool <-
function(sampleId = 1:length(area), area, targetArea, maxTry = 100L) {
    
    # check arguments
    nSamples <- length(sampleId)
    stopifnot(length(area) == nSamples)
    
    # make sure that targetArea is an ordered vector of length 2
    targetArea <- range(targetArea)
    
    # initialize pool identifier
    poolId <- rep(x = NA_integer_, length = nSamples)
    
    # initialize sample index
    sampleIndex <- 1:nSamples
    
    # exclude samples with an area greater than the maximum target pool area
    bigArea <- area > targetArea[2]
    if (any(bigArea)) {
        sampleIndex <-  setdiff(sampleIndex, which(bigArea))
    }
    
    # check trivial cases
    if (all(area %inInterval% targetArea)) {
        message("No pooling necessary")
        return(sampleIndex)
    }
    if (all(area > targetArea[2])) {
        stop("No pooling possible", call. = FALSE)
        return(NULL)
    }
    s <- sum(area)
    if (s < targetArea[1]) {
        stop("No pooling possible", call. = FALSE)
        return(NULL)
    }
    if (s %inInterval% targetArea) {
        poolId[] <- 1L
        return(poolId)
    }
    

    # estimate the (theoretical) maximum number of samples in a pool
    # Note: the theoretical minimum is not necessary as cumsum is 
    # also computed in the repeat-loop below.
    pooledArea <- cumsum(sort(area))
    index <- which(pooledArea > targetArea[2])
    maxSamples <- min(index)
                  
    # initialize number of pools
    nPools <- 0L
    
    # start pooling...
    nTry <- maxTry
    repeat {

        # select samples _without_ replacement 
        # (safer than 'sample'-function: no undesired behavior when n = 1)
        index <- sampleIndex[
            sample.int(
                n = nSamples, 
                size = min(nSamples, maxSamples), 
                replace = FALSE
            )
        ]

        # selected smallest set that is in the specified range of target areas
        pooledArea <- cumsum(area[index])
        inRange <- pooledArea %inInterval% targetArea
        if (any(inRange)) {

            # reset nTry
            nTry <- maxTry

            poolSize <- min(which(inRange))
            index <- index[1:poolSize]
            
            # construct new pool identifier
            nPools <- nPools + 1L
            poolId[index] <- nPools
            
            # update sample indices accordingly
            sampleIndex <- setdiff(sampleIndex, index)
            
            # update remaining number of samples to pool
            nSamples <- nSamples - poolSize
            
        } else {
            # termination criterion:
            # - remaining number of samples is smaller than maximum samples in a pool
            # - remaining sample area is smaller than minimum target area.
            if (nSamples <= maxSamples) {
                if (pooledArea[length(pooledArea)] < targetArea[1]) {
                    break
                }
            }
            
            # termination criterion:
            # decreases and test nTry (trial-and-error)
            nTry <- nTry - 1L
            if (nTry == 0L) {
                # give up
                break
            }
        }
        
        # termination criterion:
        if (nSamples == 0L) {
            break
        }
    }
    
    # return pool identifiers
    poolId
}




#'  Ecological Quality Ratio (EQR)
#'
#'  The ecological quality ratio is the ratio beween a parameter value and its
#'  reference value:  
#'  \deqn{EQR = \frac{x-bad}{good-bad}}{EQR = (x-bad)/(good-bad)}
#'  Depending on \code{bad} and \code{good}, the EQR usually 
#'  (but not necessarily!) varies between 0 (bad ecological quality) and 1 (good
#'  ecological quality).
#'  
#'
#' 	@param x numeric vector containing benthic indices
#'  @param bad the reference value for a bad status
#'  @param good the reference value for a good status
#'
#'  @return  numeric vector with EQR-values: low values indicate bad ecological
#'      quality and high values indicate good ecological quality.
#'  
#'  @export
eqr <- function(x, bad, good) {
    stopifnot((length(bad)  == 1L) | length(bad)  == length(x))
    stopifnot((length(good) == 1L) | length(good) == length(x))
    if (all(good > bad) | all(good < bad)) {
        return((x - bad) / (good - bad))
    }
    stop(
        paste0(
            "reference values for 'good' and 'bad' status are inconsistent.\n",
            "'good' should always be either smaller or greater than 'bad'."
        ), 
        call. = FALSE
    )
}



#'  Renaming Taxon Names
#'
#'  Convert taxon name \code{x} to taxon name \code{to} by looking it up in
#'  \code{from}. Look-up is case insensitive.
#'
#'  @param x character vector with names
#'  @param from character vector of old names
#'  @param to character vector of new names
#'  
#'
#'  @return character vector of \code{length(x)} with converted names
rename <-
function(x, from, to) {
    
    # check arguments
    if (length(from) != length(to)) {
        stop("'from' should have the same length as 'to'.", call. = FALSE)
    }

    # strip 'sp.' from TAXON name
    x <- sub(pattern = " sp\\.$", replacement = "", x = x, ignore.case = TRUE)
    
    # initialize result
    r <- rep.int(x = NA_character_, times = length(x))
    
    # keep 'x' when it is already in 'to'
    r <- to[match(x = tolower(x), table = tolower(to))]
    
    # convert 'x' to 'to' if it is in 'from'
    index <- which(is.na(r))
    if (length(index) > 0L) {
        r[index] <- to[match(x = tolower(x[index]), table = tolower(from))]
    }
    
    # return result
    r
}


#'  Harmonize Case
#'
#'  Convert text to the most occuring case. In case of ties, the first
#'  occurence in sorted order will be taken.
#'
#'  @param x character vector
#'
#'  @return character vector with harmonized names (i.e., same case)
#'  
#'  @examples 
#'  x <- c("FOO", "Foo", "bar", "FOO", "bar", "FOO", "Bar")
#'  y <- BEQI2:::harmonize(x)
#'  stopifnot(all.equal(y, c("FOO", "FOO", "bar", "FOO", "bar", "FOO", "bar")))
#'  
harmonize <- function(x) {
    lx <- tolower(x)
    lut <- tapply(
        X = x, 
        INDEX = lx, 
        FUN = function(x) {
            f <- tapply(X = x, INDEX = x, FUN = length)
            names(f)[which.max(f)]
        }
    )
    as.character(lut[match(x = lx, table = names(lut))])
}



#'  Abundance
#'
#' 	Compute abundance for each taxon.
#'
#' 	@param taxon character vector with taxa
#'  @param count integer vector with counts
#'  @export
#'
#'  @return integer vector with abundance per taxon.
abundance <- 
function(taxon, count) {
    
    # check arguments
    stopifnot(length(taxon) == length(count))

    # remove azoic (sub)samples as they do not affect abundance and
    # should be removed for further processing 
    # (e.g. species richness) anyway.
    index <- which(isAzoic(x = taxon))
    if (length(index) > 0L) {
        taxon <- taxon[-index]
        count <- count[-index]
        if (length(count) == 0L) {
            return(0L)
        }
    }

    # compute abundance
    tapply(X = count, INDEX = taxon, FUN = sum)
}



#'  Species Richness 
#'
#'  Species richness (\eqn{S}{S}) is defined as the number of taxa 
#'  (lowest identification level possible) per sampling unit 
#'  (data pool or box core sample).
#'
#' 	@param taxon character vector with taxa
#'  @param count integer vector with counts
#'  @export
#'
#'  @return species richness (integer vector of length 1)
speciesRichness <- 
function(taxon, count) {

    # abundance (delegate argument checking and removal of azoic samples)
    n <- abundance(taxon = taxon, count = count)
    
    # species richness (after removal of azoic samples)
    length(n)
}


#'  Margalef Index of Diversity
#'
#'  Margalef Index of Diversity is given by
#'  \deqn{D = \frac{S-1}{\ln(N)}}{D = (S-1)/ln(N)}
#'  
#'  For \eqn{N=1}{N=1}, the index is set to 0.
#'
#'  @param taxon character vector with taxa
#'  @param count integer vector with counts
#'  @export
#'
#'  @return Margalef diversity index (numeric vector of length 1)
margalef <- 
function(taxon, count) {

    # species richness
    S <- speciesRichness(taxon = taxon, count = count)

    # abundance (delegate argument checking and removal of azoic samples)
    n <- abundance(taxon = taxon, count = count)
    
    # total abundance (azoic samples are excluded)
    N <- sum(n)
    
    # Margalef's index of diversity (azoic samples are excluded)
    if (N == 1L) {
        D <- 0 
    } else {
        D <-(S - 1) / log(N)
    }
    
    # return D
    D
}


#'  Shannon's Entropy 
#'
#'  Compute entropy according to Shannon (1948)
#'
#'  @param taxon character vector with taxa
#'  @param count integer vector with counts
#'  @export
#'
#'  @return Shannon's entropy
#'  
#'  @references Shannon, C. E., 1948. A Mathematical Theory of Communication.
#'      Bell System Technical Journal 27: 379-423.
entropy <- 
function(taxon, count) {
    n <- abundance(taxon = taxon, count = count)
    p <- n / sum(n)
    -sum(p * log2(p))
}



#'  Create BEQI-2 Directory Structure
#'
#'  Creates a BEQI2-directory structure and populates it with some
#'  relevant BEQI2-files. Users may wish to modify this directory structure 
#'  and add their own data.
#'
#'  @param path name of an exisiting directory. This directory should
#'      be empty to prevent loss of data. If missing, a dialogue will
#'      appear.
#'  
#'  @export
BEQI2dir <- function(path = NULL) {

    # interactive selection
    if (is.null(path)) {
        if (capabilities("tcltk")) {
            path <- tk_choose.dir(
                caption = "Select directory to store BEQI-2 files"
            )
        } else {
            stop(
                "The 'tcltk'-package is not supported on this machine.\n",
                "Please provide a valid path as function argument\n",
                call. = FALSE
            )
        }
    }
    
    # check path
    if (!file.exists(path)) {
        stop("directory does not exist", call. = FALSE)
    }
    
    # check if directory is empty (to prevent overwriting existing data)
    if (length(list.files(path)) != 0L) {
        stop(sprintf("directory %s is not empty!", sQuote(path)), call. = FALSE)
    }
    
    # populate directories
    tmp <- file.copy(
        from = c(
            system.file("extdata/INPUT-FILES", package = "BEQI2"),
            system.file("extdata/REF-FILES", package = "BEQI2"),
            system.file("extdata/settings.json", package = "BEQI2")
        ),
        to = path,
        recursive = TRUE
    )

    # show message
    message(
        sprintf(
            paste0(
                "Directory %s\nhas been populated with BEQI-2 files.\n",
                "To run the BEQI-2 tool, type: BEQI2()\n",
                'For the tutorial, type: vignette("BEQI2")\n',
                "For more technical information, type: ?beqi2"
            ), 
            sQuote(path)
        )
    )
}



#'  Genus to Species Conversion
#'
#' 	For each sample, the algorithm tries to convert taxa on the genus level
#'  to the species level. Counts at the genus level will be distributed 
#'  over the species level proportional to the available number of species.
#'  
#'  @note The updated counts are not necessarily integers.
#'
#' 	@param id sample identifier
#'  @param taxon taxon name on either the genus or species level
#'  @param count total number of individuals of a specific taxon in a sample
#'
#'  @return \code{data.frame} with three columns: \code{id} the sample
#'      identifier, \code{taxon} the taxon name, and \code{count} the  count
#'      after genus to species conversion.
#'  
#'  @export
#'      
#'  @examples
#'      genusToSpecies(id = c(1, 1, 1), 
#'          taxon = c("Genus1", "Genus1 s1", "Genus1 s2"), count = c(4, 2, 1))
genusToSpecies <-
function(id, taxon, count) {
    
    # check arguments
    stopifnot(length(id) == length(taxon))
    stopifnot(length(id) == length(count))
    taxon <- as.character(taxon)

    # split taxon into a genus and a species part
    d <- strsplit(x = taxon, split =" +")
    genusPart <- sapply(X = d, FUN = function(x){x[1]}) 
    speciesPart <- sapply(X = d, FUN = function(x) {
        if (length(x) == 1L) {
                NA_character_
            } else {
                paste(x[-1], collapse = " ")
            }
        }
    )
    
    # convert genus to proportionally to existing species
    d <- ddply(
        .data = data.frame(id, taxon, genusPart, speciesPart, count), 
        .variables = c("id", "genusPart"), 
        .fun = function(x) {
            isGenusLevel <- is.na(x$speciesPart)
            isSpeciesLevel <- !isGenusLevel
            if (any(isGenusLevel) & any(isSpeciesLevel)) {

                # extract counts of species
                n <- x$count[isSpeciesLevel]
                
                # convert counts to weights
                w <- n / sum(n)
                
                # assign genus counts proportionally to species
                p <- w * x$count[isGenusLevel]
    
                # update counts
                x$count[isGenusLevel] <- 0L
                x$count[isSpeciesLevel] <- n + p
            }
            x
        }
    )

    # return result
    d[, c("id", "taxon", "count")]
    
}



