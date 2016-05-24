guiSuiteList <- function (sep = "\t", path = NULL, compare = TRUE)
{
    Suites <- svSuiteList()
	if (compare) {
		oldSuites <- .getTemp(".guiSuiteListCache", default = "")
		## Compare both versions
		if (!identical(Suites, oldSuites)) {
			## Keep a copy of the last version in SciViews:TempEnv
			.assignTemp(".guiSuiteListCache", Suites)
			Changed <- TRUE
		} else Changed <- FALSE
	} else {
		Changed <- TRUE
		.assignTemp(".guiSuiteListCache", Suites)
	}
    if (is.null(path)) { # Return result, as a single character string with sep
		if (Changed) {
			if (!is.null(sep)) Suites <- paste(Suites, collapse = sep)
			return(Suites)
		} else return(NULL)
	} else { # Write to a file called 'Suites.txt' in this path
		file <- file.path(path, "Suites.txt")
		if (Changed) {
			if (is.null(sep)) sep <- "\n"
			cat(Suites, sep = sep, file = file)
		}
		return(invisible(Changed))
	}
}

guiSuiteAutoList <- function (...)
{
	## Is koCmd() available?
	if (!exists("koCmd", mode = "function")) return(TRUE)
	## Is it something changed in the unit list?
	res <- guiSuiteList(sep = ",", path = NULL, compare = TRUE)
	if (!is.null(res))
		ret <- get("koCmd")('sv.r.unit.getRUnitList_Callback("<<<data>>>");',
			data = res)
	return(TRUE)
}

guiTestFeedback <- function (object, path = NULL, ...)
{
	## Give feedback to client about the currently running tests
	## TODO: feedback about test run
}

guiTestReport <- function (object, sep = "\t", path = NULL, ...)
{
	## Report the results of a test to the GUI client
	if (!is.svSuiteData(object))
		stop("'object' must be a 'svSuiteData' object")

	## For all 'svTestData' objects, create a table with test results for the GUI
	## Indicate global results of the Unit Test
	Tests <- ls(object)
    if (length(Tests) == 0) {
        Res <- "<<<svUnitSummary>>>|||0|||0|||0|||0"
    } else {
        ## Get general information about the tests
        Stats <- stats(object)
		Tests <- rownames(Stats)	# To make sure we use the same!
		Stats$label <- paste(">", sub("^test", "", Tests), " (",
			round(Stats$timing, 3), " sec)", sep = "")
		State <- table(Stats$kind)
		Res <- paste("<<<svUnitSummary>>>|||", State[1], "|||", State[2],
			"|||", State[3], "|||", State[4], sep = "")
		Kinds <- as.numeric(Stats$kind)
		Kinds[Kinds == 4] <- 0 	# Use 0 instead of 4 for deactivated tests
        Stats$kind <- Kinds
		## Get the type for the objects
		Units <- Stats$unit
		Types <- rep("units in packages", length(Units))
		Types[Units == ""] <- "other objects"
		## TODO: include also dirs!
		Dir1 <- gsub("\\\\", "/", dirname(Units))
		Dir2 <- dirname(Dir1)
		Dir3 <- dirname(Dir2)
		TempDir <- gsub("\\\\", "/", tempdir())
		Types[Dir1 == TempDir] <- "objects in .GlobalEnv"
		Types[tolower(basename(Dir2)) == "inst" ||
			tolower(basename(Dir3)) == "inst"] <- "units in sources"
		## Keep only "*" in Units
		Units <- basename(Units)
		Units[regexpr("^runit.+\\.[rR]$", Units) == -1] <- ""
		Units[Dir1 == TempDir] <- "" # No second level for objects in .GlobalEnv
		Units <- sub("^runit(.+)\\.[rR]$", "\\1", Units)
		change <- Units != ""
		Units[change] <- paste(">unit", Units[change])
		## Complete label is Type<Unit<Test (x.xxx sec)
		Stats$label <- paste(Types, Units, Stats$label, sep = "")
		## Sort Tests and Stats according to label alphabetically
		ord <- order(Stats$label)
		Stats <- Stats[ord, ]
		Tests <- Tests[ord]
		## Get detailed information about each test
        lastUnit <- ""
		for (Test in Tests) {
			Data <- Stats[Test, ]
			## Calculate Info
			tData <- Log()[[Test]]
			tStats <- stats(tData)
			Info <- paste(c("Pass:", "Fail:", "Errors:"), tStats$kind[1:3],
            collapse = " ")
			## Don't print tests that succeed if !all
			tData <- tData[tData$kind != "OK", ]
			## Get info about each individual filtered test
			if (nrow(tData) > 0) {
				Result <- ifelse(tData$res == "", "",
					paste("\n", tData$res, sep = ""))
				Info <- paste(Info, "\n", paste("* ", tData$msg, ": ",
					tData$call, .formatTime(tData$timing, secDigits = 3),
					" ... ", as.character(tData$kind), Result, sep = "",
					collapse = "\n"), sep = "")
			}
			## Calculate URI (currently, the name of the unit file
			## and the name of the test function)
			if (Data$unit == "") URI <- Data$unit else
				URI <- paste(Data$unit, Test, sep = "#")
			if (Data$unit != lastUnit) {
				lastUnit <- Data$unit
				Res <- c(Res, paste("<<<svUnitFile>>>|||", Data$unit,
					"|||||||||", sep = ""))
			}
			Res <- c(Res, paste("<<<svUnitTest>>>|||", Data$label, "|||",
				Data$kind, "|||", Info, "|||", URI, sep = ""))
		}
	}
	Res <- paste(gsub("\t", "    ", Res), collapse = sep)
	if (is.null(path)) {
		return(Res)
	} else {
		cat(Res, file = path)
	}
    return(path)
}
