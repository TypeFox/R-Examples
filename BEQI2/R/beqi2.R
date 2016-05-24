#' 	Perform BEQI-2 Analysis
#'
#' 	This function performs a complete BEQI-2 analysis following the
#'	settings provided in \code{filename}.
#'
#' 	@param filename name of the JSON file defining all analysis steps.
#'	@param tmpdir directory to store temporary files (for debugging only)
#'  @param browse load resulting report in a browser? \code{TRUE} or \code{FALSE}
#'
#' 	@export
beqi2 <-
function(filename = NULL, tmpdir = tempfile(pattern = "BEQI2"), browse = TRUE) {

    # prevent potential problems with dates in other locales
    oldLocale <- Sys.getlocale("LC_TIME")
    on.exit(Sys.setlocale("LC_TIME", oldLocale))
    Sys.setlocale("LC_TIME", "C")
    
    # interactive selection of filename
    if (is.null(filename)) {
        if (capabilities("tcltk")) {
            filename <- tk_choose.files(
                default = "", 
                caption = "Select file with BEQI2 settings",
                multi = FALSE, 
                filters = matrix(data = c("BEQI2 settings", ".json"), nrow = 1)
            )
        } else {
            stop(
                "The 'tcltk'-package is not supported on this machine.\n",
                "Please provide a valid filename as function argument\n",
                call. = FALSE
            )
        }
    }

    # stop if the user presses Cancel or Esc
    if(length(filename) == 0L) {
        message("The BEQI-2 run has been cancelled by the user.")
        return(invisible(NULL))
    }
    
    # check if filename exists
    if (!file.exists(filename)) {
        stop(
            sprintf("JSON-file %s does not exist", sQuote(filename)), 
            call. = FALSE
        )
    }

    # initialization message
    message("The BEQI-2 tool is running...")

    # read settings
	settings <- readSettings(filename)
    
    # set working directory
    owd <- getwd()
    on.exit(setwd(owd))
    setwd(dirname(filename))

    # normalize paths (full paths to make package more robust)
    for (f in names(settings$files)) {
        settings$files[[f]] <- suppressWarnings(normalizePath(settings$files[[f]]))
    }
    
    # add output files
    outputDir <- file.path(
        getwd(), 
        paste0("OUTPUT-", format(Sys.time(), format = "%Y%m%dT%H%M%S"))
    )
    dir.create(outputDir)
    prefix <- sub(
        pattern = "\\.[[:alnum:]]+$", 
        replacement = "", 
        x = basename(settings$files$beqi2)
    )
    settings$files$log <- file.path(outputDir, 
        paste0("LOG-", prefix, ".log"))
    settings$files$out_ecotope <- file.path(outputDir, 
        paste0("ECOTOPE-", prefix, ".csv"))
    settings$files$out_objectid <- file.path(outputDir, 
        paste0("OBJECTID-", prefix, ".csv"))
    if (settings$pooling$enabled) {
        settings$files$pooling <- file.path(outputDir, 
            paste0("POOLING-", prefix, ".csv"))
    }
    settings$files$report <- file.path(outputDir, 
        paste0("REPORT-", prefix, ".html"))

    # start log-file
    toLog <- function(level = c("INFO", "WARNING", "ERROR"), message) {
        level <- match.arg(level)
        cat(
            format(Sys.time()), " [", level, "] ", message, "\n", 
            sep = "",
            file = settings$files$log, 
            append = TRUE
        )
        if (level != "INFO") {
            message <- paste(message, "(see log-file)")
            switch(level,
               "ERROR"   = stop(message, call. = FALSE),
               "WARNING" = warning(message, call. = FALSE)
            )
        }
    }
    toLog("INFO", "Starting a new BEQI-2 session")
    on.exit(toLog("INFO", "This BEQI-2 session has been terminated"), add = TRUE)

    # initialize random number generator
    if (settings$pooling$enabled) {
        toLog("INFO", "Initializing the pseudo random number generator...")
        if (is.null(settings$pooling$randomseed)) {
            toLog("INFO", "No seed has been specified.")
            toLog("INFO", "The default initialization process will be followed.")
        } else {
            set.seed(seed = settings$pooling$randomseed)
        }
        toLog("INFO", "the pseudo random number generator has been initialized.")
    }
    
    
    # check existence of BEQI-2 file
    toLog("INFO", 
          sprintf(
                "Checking the existence of BEQI-2 file %s...",
                sQuote(basename(settings$files$beqi2))
          )
    )
    if (!file.exists(settings$files$beqi2)) {
        toLog("ERROR", "The BEQI-2 file has not been found")
        return(invisible(NULL))
    }
    toLog("INFO", "the BEQI-2 file has been found")

    # read BEQI2-file
    toLog("INFO", "Reading the BEQI2-file...")
    d_beqi <- tryCatch(
        readBEQI(filename = settings$files$beqi2),
        error = function(e) {
            toLog("ERROR", sprintf("while reading BEQI2-file. %s", e$message))
        }
    )
    toLog("INFO", "the BEQI-2 file has been read")
    
    # check if records are within the period of interest
    month <- as.integer(format(d_beqi$DATE, format = "%m"))
    inPOI <- (month >= settings$months[1]) & 
             (month <= settings$months[2])
    if (!any(inPOI)) {
        toLog("ERROR", 
            sprintf("No months in file %s are in the specified interval [%s].",
                sQuote(basename(settings$files$beqi2)),
                paste(settings$months, collapse = ", ")
            )
        )
    }

    # check existence of TWN-file
    toLog("INFO", 
          sprintf(
                "Checking the existence of TWN-file %s...",
                sQuote(basename(settings$files$speciesnames))
          )
    )
    if (!file.exists(settings$files$speciesnames)) {
        toLog("ERROR", "The TWN-file has not been found")
        return(invisible(NULL))
    }
    toLog("INFO", "the TWN-file has been found")

    # read TWN-file
    toLog("INFO", "Reading the TWN-file...")
    d_twn <- tryCatch(
        readTWN(filename = settings$files$speciesnames),
        error = function(e) {
            toLog("ERROR", sprintf("while reading TWN-file. %s", e$message))
        }
    )
    toLog("INFO", "the TWN-file has been read")
    
    # read AMBI's
    d_sens <- suppressMessages(readAMBI())
    names(d_sens) <- c("TAXON", "AMBI")
    d_sens$TAXON <- rename(x = d_sens$TAXON, from = d_twn$taxonname, to = d_twn$taxon)
    d_sens <- unique(na.omit(d_sens))
    if (!is.null(settings$files$ambi) && (settings$files$ambi != "")) {
        toLog("INFO", 
              sprintf(
                    "Checking the existence of AMBI-file %s...",
                    sQuote(basename(settings$files$ambi))
              )
        )
        if (!file.exists(settings$files$ambi)) {
            toLog("ERROR", "the AMBI-file has not been found")
            return(invisible(NULL))
        }
        toLog("INFO", "the AMBI-file has been found")
        toLog("INFO", "Reading the AMBI-file...")
        d <- tryCatch(
            readAMBI(filename = settings$files$ambi),
            error = function(e) {
                toLog("ERROR", sprintf("while reading AMBI-file. %s", e$message))
            }
        )
        toLog("INFO", "the AMBI-file has been read")
        names(d) <- c("TAXON", "AMBI_user")
        d$TAXON <- rename(x = d$TAXON, from = d_twn$taxonname, to = d_twn$taxon)
        d <- unique(na.omit(d))
        d <- merge(x = d_sens, y = d, all = TRUE) 
        sel <- is.na(d$AMBI_user)
        d$AMBI_user[sel] <- d$AMBI[sel]
        d_sens <- data.frame(TAXON = d$TAXON, AMBI = d$AMBI_user,
                             stringsAsFactors = FALSE)
    }

    # read ecotope reference file
    toLog("INFO", 
          sprintf(
                "Checking the existence of ecotope reference file %s...",
                sQuote(basename(settings$files$ecotopes))
          )
    )
    if (!file.exists(settings$files$ecotopes)) {
        toLog("ERROR", "the ecotope reference file has not been found")
        return(invisible(NULL))
    }
    toLog("INFO", "the ecotope reference file has been found")
    toLog("INFO", "Reading the ecotope reference file...")
    extra <- NULL
    if (!is.null(settings$files$iti)) {
        extra <- c(extra, "ITI")
    }
    if (!is.null(settings$files$fibi)) {
        extra <- c(extra, "FIBI")
    }
    d_erf <- tryCatch(
        readERF(filename = settings$files$ecotopes, extra = extra),
        error = function(e) {
            toLog("ERROR", sprintf("while reading ecotope reference file. %s", e$message))
        }
    )
    toLog("INFO", "the ecotope reference file has been read")

    # check if reference data are available for all records in d_beqi
    toLog("INFO", "Checking if reference data are available for all records in the BEQI2-file...")
    tmp1 <- unique(paste(d_beqi$OBJECTID, d_beqi$ECOTOPE, sep = "/"))
    tmp2 <- unique(paste( d_erf$OBJECTID,  d_erf$ECOTOPE, sep = "/"))
    isMissing <- is.na(match(x = tmp1, table = tmp2))
    if (any(isMissing)) {
        toLog("ERROR", sprintf(
                "The following combinations of %s are missing in the ecotope reference file: %s",
                sQuote("OBJECTID/ECOTOPE"),
                toString(tmp1[isMissing])
            )
        )
        return(invisible(NULL))
    }
    toLog("INFO", "reference data are available for all records in the BEQI2-file.")
    
	# create temporary directory
	if (!file.exists(tmpdir)) {
        toLog("INFO", "Creating a temporary directory...")
		dir.create(tmpdir)
	}
    toLog("INFO", "a temporary directory has been created.")

	# copy template of the report to temporary directory
    toLog("INFO", "Populating the temporary directory...")
	templates <- list.files(
		path = system.file("Rmd", package = "BEQI2"),
		pattern = "\\.Rmd$", full.names = TRUE)
	file.copy(from = templates, to = tmpdir)
    toLog("INFO", "the temporary directory has been populated.")

    # create Markdown document 
    # (code below works better than knit2html)
    toLog("INFO", "Starting to create a report...")
    setwd(tmpdir)
	suppressMessages(
        mdfile <- try(knit(input = "beqi2.Rmd", quiet = TRUE), silent = TRUE)
    )
    if (inherits(mdfile, "try-error")) {
        toLog(
            level = "ERROR", 
            message = toString(attr(mdfile, "condition")$message)
        )
        return(invisible(NULL))
    }
    toLog("INFO", "a report has been created.")
    toLog("INFO", "Converting the report to HTML...")
	output <- markdownToHTML(
		file  = mdfile, 
		output = NULL,
        options = getOption("markdown.HTML.options"),
        extensions = getOption("markdown.extensions"),
    	title = "Benthic Ecosystem Quality Index 2 Report",
        stylesheet = system.file("css", "beqi2.css", package = "BEQI2")
	)
    writeLines(text = output, con = settings$files$report)
    toLog("INFO", "the report has been converted to HTML.")

	# view result
	if (browse) {
		browseURL(settings$files$report)
	}
    
    # finalization
    message("The BEQI-2 run has been completed successfully.")
}



#'  Perform BEQI-2 Analysis
#'
#'  @inheritParams beqi2
#'  
#'  @rdname beqi2
#'  
#' 	@export
BEQI2 <- 
function(filename = NULL, tmpdir = tempdir(), browse = TRUE) {
    beqi2(filename = filename, tmpdir = tmpdir, browse = browse)
}
