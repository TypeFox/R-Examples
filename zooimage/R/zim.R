## Copyright (c) 2004-2012, Ph. Grosjean <phgrosjean@sciviews.org>
##
## This file is part of ZooImage
##
## ZooImage is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## ZooImage is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with ZooImage.  If not, see <http://www.gnu.org/licenses/>.

## Functions for manipulating .zim files (ZooImage Metadata/measurements)
## These.zim files contain metadata required to analyze plankton images
## and to record the way they were processed. Measurements on each identified
## object can also be appended in a table at the end of this file (in this case,
## the usual extension is '_dat1.zim' to indicate that data processed with
## ZooImage version 1 are present in the file).

## Check if a file is a "(_dat1).zim" file (must start with "ZI1", "ZI2"
## or "ZI3" and have a '.zim' extension)
isZim <- function (zimfile)
{
	## Check if the file does not exist or is a directory
	if (!checkFileExists(zimfile, force.file = TRUE, extension = "zim"))
		return(FALSE)

	## Check the first line
	return(checkFirstLine(zimfile))
}

## Make required .zim files for one or more images
zimMake <- function (dir = ".", pattern = extensionPattern("tif"),
images = list.files(dir, pattern))
{
	## Check that there are images to process
	if (!length(images)) {
		warning("no images to process!")
		return(invisible(FALSE))
	}

	## Name of images is something like SCS.xxxx-xx-xx.SS+Ann.tif
	## We make the same .zim file for all ...+Ann images, so, reduce the list
	zims <- sort(unique(sampleInfo(images, type = "fraction",
		ext = pattern)))
	zims <- file.path(dir, paste0(zims, ".zim"))
		
	## The process to run in batch
	makeZim <- function (zim) {
		if (!file.exists(zim)) {
			message("Processing ZIM file ", basename(zim), "...")
			zimCreate(zim, template = getTemp(".template"), wait = TRUE)
			## Get this zim as new template
			assignTemp(".template", zim)
		} else message("Checking ZIM file ", basename(zim), "...")
		## Verify that the .zim file is correct
		zimVerify(zim) >= 0
	}
	on.exit(rmTemp(".template"))
	message("Creating ZIM files...")
	flush.console()
	ok <- batch(zims, makeZim, verbose = FALSE)
	if (!ok) {
		warning(sum(attr(ok, "ok")), "/", length(images),
			" ZIM files created (see .last.batch)")
		invisible(FALSE)
	} else {
		message("-- Done! --")
		invisible(TRUE)
	}
}

## Verify a "(_dat1).zim" file (all required fields + return the number of items
## in it). If it succeeds, return the number of measured items as numeric
## Otherwise, return -1... If there is no data, return 0
zimVerify <- function (zimfile, is.dat1 = hasExtension(zimfile, "_dat1.zim"),
check.table = FALSE)
{
	## Required fields
	## Here are predefined required fields before measurements
	reqfields <- c("[Image]", "Author", "Hardware", "Software",
        "ImageType", "[Fraction]", "Code", "Min", "Max", "[Subsample]",
        "SubPart", "SubMethod", "CellPart", "Replicates", "VolIni",
        "VolPrec")

	## Then required fields when measurements are done
	reqfields2 <- c("[Process]")
	## Finally, required column headers
    reqcols <- c("!Item", "Label", "BX", "BY", "Width", "Height")

	## Determine if there are custom verification rules defined and if they
	## are active
    newRules <- getOption("ZI.zim")
    if (length(newRules) && newRules$active == TRUE) {
        ## Do we delegate the whole process to a custom verification function?
		verifyAll <- newRules$verify.all
        if (!is.null(verifyAll) && inherits(verifyAll, "function"))
            return(verifyAll(zimfile = zimfile, is.dat1 = is.dat1,
				check.table = check.table))

		## Should we use additional verification code instead?
		verify <- newRules$verify
        reqfields <- c(reqfields, newRules$zim.required)
        reqfields2 <- c(reqfields2, newRules$dat1.zim.required)
        reqcols <- c(reqcols, newRules$dat1.data.required)
    } else verify <- NULL

	## Check that it is a zimfile
	if (!isZim(zimfile)) return(-1)

	## Run first the extra verification code
	if (!is.null(verify) && inherits(verify, "function")) {
		## We need to grab the error here and call stop from here to maintain
		## the old API and to allow the custom version of stop to be called
		## with the correct context of the zimVerify() function
		res <- try(verify(zimfile, is.dat1 = is.dat1,
            check.table = check.table), silent = TRUE)
		if (inherits(res, "try-error") || (is.character(res) && nchar(res))) {
			warning(as.character(res))
			return(-1)
		}
    }

	## Read the file...
	## Equal sign is used as comment char, in order to read only the field names
    Lines <- scan(zimfile, character(), sep = "\t", skip = 1, flush = TRUE,
		quiet = TRUE, blank.lines.skip = FALSE, comment.char = "=")
	## Trim leading and trailing white spaces
	Lines <- trimString(Lines)

	## Check that all required fields are present for a simple .zim file
    misfields <- reqfields[!(reqfields %in% Lines)]
    if (length(misfields) > 0) {
        warning(paste("Missing fields:", paste(misfields, collapse = ", ")))
		return(-1)
	}

	## Check if this is a _dat1.zim file with measurements
    if ("[Data]" %in% Lines) {
        ## Check for missing fields
		misfields2 <- reqfields2[!(reqfields2 %in% Lines)]
        if (length(misfields2) > 0) {
            warning(paste("Missing [Process] fields:", paste(misfields2,
				collapse = ", ")))
			return(-1)
		}

		## Check for required column headers
		posHeaders <- grep("^\\[Data\\]$", Lines)[1] + 1
		LineHeader <- scan(zimfile, character(), sep = "%", skip = posHeaders,
			nmax = 1, flush = TRUE, quiet = TRUE, comment.char = "=")
		Headers <- trimString(strsplit(LineHeader, "\t")[[1]])
		misHeaders <- reqcols[!(reqcols %in% Headers)]
		if (length(misHeaders) > 0) {
		    warning(paste("Missing columns in the table:", paste(misHeaders,
				collapse = ", ")))
			return(-1)
		}

		## Check that the table can be read
        if (isTRUE(as.logical(check.table))) {
			## Check the [Data] section
            posMes <- grep("^\\[Data\\]$", Lines)
            if (length(posMes) == 0) {
                warning("Trying to read the table of measurements but no [Data] section found!")
				return(-1)
            } else { # The [Data] section is found
				## we try to call read.table, catch the error, and throw it again
				## from here, because stop might have a different meaning
				## in the context of the zimVerify() function
				## allowing to use the zooImage calling handlers,
				## see errorHandling.R
				Mes <- try(read.table(zimfile, sep = "\t", header = TRUE,
					skip = posMes + 1), silent = TRUE)
				if (inherits(Mes, "try-error")) {
					warning(paste("Unable to read the table of measurements! : ",
						Mes))
					return(-1)
				} else { 	# Successful reading of the table of measurements
					return(nrow(Mes))	# Return the number of items measured
				}
            }
        } else {
			## Alternative method that does not read the table
			## We don't read the table, use a different method to get the number
			## of entries in it
			## Read the last entry in Lines and convert it to a numeric value:
			## should be the number of items measured
			nItems <- Lines[length(Lines)]
			if (sub("^[0-9]+$", "", nItems) != "") {
			    warning("Impossible to determine the number of items measured!")
				return(-1)
			}
			return(as.integer(nItems))
        }
    } else if (isTRUE(as.logical(is.dat1))) {
		warning("No measurements found in this file")
		return(-1)
	} else return(0)
}

## Extract notes from .zip files and place them in .zim files
zimExtractAll <- function (zipdir = ".", zipfiles = zipList(zipdir),
path = NULL, replace = FALSE)
{
	## Make sure all zipfiles are in the same directory
	zipdirs <- dirname(zipfiles)
	if (length(unique(zipdirs)) > 1) {
		warning("all ZIP files must be located in the same directory!")
		return(invisible(FALSE))
	}

	## Check that the dir exists!
	if (!checkDirExists(zipdir)) return(invisible(FALSE))

	## Move to zipdir
	initdir <- setwd(zipdir)
	on.exit(setwd(initdir))
	zipdir <- getwd()   # That way, if we had ".", it is now expanded

	## Use only basenames for zip files
	zipfiles <- sort(basename(zipfiles))

	## Check that zipfiles exist
	if (!all(file.exists(zipfiles))) {
		stop("one or several ZIP files not found!")
		return(invisible(FALSE))
	}

	## Look at the path where to place .zim files
	if (!length(path)) {
		## The rule is the following one:
		## 1) if last subdir is "_raw", then place .zim file up one level
		## 2) else, place them in the same dir as the zip files
		path <- zipdir
		if (tolower(basename(path)) == "_raw") path <- dirname(path)
	} else {    # Look if this is a valid directory
		path <- path[1]
		if (!checkDirExists(path)) return(invisible(FALSE))
	}

	## Compute the names of .zim files from the names of .zip files
	## Note: use only the fraction, that is, SCS.xxxx-xx-xx.SS+F from
	## SCS.xxxx-xx-xx.SS+Fnn)
	## If there are duplicates, only extract first one
	zimfiles <- paste(sampleInfo(zipfiles, "fraction",
		ext = extensionPattern(".zip")), "zim", sep = ".")
	keep <- !duplicated(zimfiles)
	zimfiles <- zimfiles[keep]
	zipfiles <- zipfiles[keep]

	## Make full path name for zimfiles
	zimfiles <- file.path(path, zimfiles)

	## If replace == FALSE, eliminate existing .zim files from the list
	if (!isTRUE(as.logical(replace))) {
		keep <- !file.exists(zimfiles)
		zimfiles <- zimfiles[keep]
		zipfiles <- zipfiles[keep]
	}

	## Are there files left
	if (!length(zimfiles)) {
		message("done: no file to process!")
		return(invisible(TRUE))
	}

	## Extract .zim files, one at a time, and check them
	items <- 1:length(zimfiles)
	zimExtract <- function (item, zipfiles, zimfiles) {
		## Extract the .zim file
		if (is.null(zipNoteGet(zipfiles[item], zimfiles[item])))
			return(FALSE)
		
		## Check that the .zim file is created and return result accordingly
		zimVerify(zimfiles[item]) >= 0
	}
	## Batch process all files
	message("Compression of images...")
	flush.console()
	ok <- batch(items, fun = zimExtract, zipfiles = zipfiles,
		zimfiles = zimfiles, verbose = FALSE)
	if (!ok) {
		warning(sum(attr(ok, "ok")), "/", length(zimfiles),
			" metadata were extracted into ZIM files (see .last.batch)")
		invisible(FALSE)
	} else {
		message("-- Done! --")
		invisible(TRUE)
	}
}

## Given a list of .zip files and a path where .zim files are located,
## update comment fields of the .zip files with latest .zim content
zimUpdateAll <- function (zipdir = ".", zipfiles = zipList(zipdir),
zimdir = NULL, check.zim = TRUE)
{
    ## Make sure we have full path for zip files
	if (zipdir == ".") zipdir <- getwd()
	zipfiles <- file.path(zipdir, zipfiles)

    ## Check that zipfiles exist
	if (!all(file.exists(zipfiles))) {
		warning("one or several ZIP files not found!")
		return(invisible(FALSE))
	}

	## Look for the path where .zim files are located
	if (!length(zimdir)) {
		## The rule is the following one:
		## 1) if last subdir of .zip files is "_raw", then .zim files
		##    should be up one level
		## 2) else, look at the same dir
		zimdir <- zipdir
		if (tolower(basename(zimdir)) == "_raw")
			zimdir <- dirname(zimdir)
	} else {    # Look if this is valid directory
		zimdir <- zimdir[1]
		if (!checkDirExists(zimdir, message = "'%s' is not a valid directory!"))
			return(invisible(FALSE))
	}

	## Switch to that dir
	initdir <- setwd(zimdir)
	on.exit(setwd(initdir))

	## Compute the names of zim files from the names of zip files
	## Note: use only the fraction, that is, SCS.xxxx-xx-xx.SS+F from
	## SCS.xxxx-xx-xx.SS+Fnn)
	## If there are duplicates, only extract first one
	zimfiles <- sprintf( "%s.zim",
		sampleInfo(zipfiles, "fraction", ext = extensionPattern("zip")))

	## Eliminate path for zimfiles
	zimfiles <- basename(zimfiles)

	## Keep only existing .zim files
	keep <- file.exists(zimfiles)
	zimfiles <- zimfiles[keep]
	zipfiles <- zipfiles[keep]

	## Are there files left?
	if (!length(zimfiles)) {
		message("done: no file to update!")
		return(invisible(TRUE))
	}

	## Check the zim files
	ok <- TRUE
	if (isTRUE(as.logical(check.zim))) {
		message("Verification of ZIM files...")
		flush.console()
		zfiles <- unique(zimfiles)
		zimCheck <- function (zim) {
			message("Verifying '", basename(zim), "' ...")
			zimVerify(zim) >= 0
		}
		ok <- batch(zfiles, zimCheck, verbose = FALSE)
	}
	if (ok) {
		message("-- Done! --")
	} else {
		warning("corrupted ZIM file(s) found, update not started (see .last.batch)")
		return(invisible(FALSE))
	}

	## If everything is OK, update comments in the zip files with the content
	## of the .zim files
	message("Update of metadata in ZIP files from ZIM data...")
	flush.console()
	items <- 1:length(zipfiles)
	updateZip <- function (item, zipfiles, zimfiles) {
		zipNoteAdd(zipfiles[item] , zimfiles[item])
	}
	ok <- batch(items, zipfiles = zipfiles, zimfiles = zimfiles,
		verbose = FALSE)
	if (!ok) {
		warning(sum(attr(ok, "ok")), "/", length(zipfiles),
			" ZIP files were updated (see .last.batch)")
		invisible(FALSE)
	} else {
		message("-- Done! --")
		invisible(TRUE)
	}
}

## Create a .zim file
zimCreate <- function (zimfile, template = NULL, edit = TRUE,
editor = getOption("fileEditor"), wait = FALSE)
{
	## Create a .zim file from a template and edit it
	if (missing(zimfile) || is.null(zimfile) || zimfile == "") {
		zimfile <- dlgInput("Give a name for the new ZIM file:",
			title = "ZIM file creation", default = "myfile.zim")$res
		if (!length(zimfile)) return(invisible(FALSE))
		if (!hasExtension(zimfile, "zim"))
			zimfile <- paste(zimfile, ".zim", sep = "")
	}

	## If the file exists, edit existing version instead
    if (file.exists(zimfile))
		if (isTRUE(as.logical(edit))) {
			return(zimEdit(zimfile, editor = editor, wait = wait))
		} else return(invisible(TRUE))

	## Look for the template
	if (!length(template))
		template <- file.path(getOption("ZITemplates"), "default.zim")
	if (!isZim(template)) return(invisible(FALSE))
	## Copy the template into the new file
	file.copy(template, zimfile)
	
	## Possibly edit this new file
	if (isTRUE(as.logical(edit))) {
		return(zimEdit(zimfile, editor = editor, wait = wait))
	} else return(invisible(TRUE))
}

## Edit a .zim file
zimEdit <- function (zimfile, editor = getOption("fileEditor"), wait = FALSE,
...)
{
	if (missing(zimfile) || !length(zimfile) || zimfile == "") {
		zimfile <- selectFile("Zim")
		if (zimfile == "") return(invisible(FALSE))
	} else if (!isZim(zimfile)) return(invisible(FALSE))
	fileEdit(zimfile, editor = editor, wait = wait, ...)
}

## Create a dat1.zim file by pooling lst and results.csv tables
zimDatMakeFlowCAM <- function (zimfile)
{
	## Check argument
	if (length(zimfile) != 1 || !is.character(zimfile)) {
		warning("you must select one ZIM file")
		return(invisible(FALSE))
	}
	if (!checkFileExists(zimfile, extension = "zim", force.file = TRUE))
		return(invisible(FALSE))

	message("compiling measurements and '", basename(zimfile), "'")
	
	## Dir containing the .zim file
	zidir <- dirname(zimfile)

	## Read list file
	## Read visual spreadsheet data (from the FlowCAM)
	visdata <- .lstRead(file.path(zidir,
		paste(basename(zidir), "lst", sep = ".")), skip = 2)
	
	## Read ImageJ results (from FITVis)
	fitdata <- read.table(file.path(zidir, "results.csv"),
		sep = ",", header = TRUE, dec = ".")
	
	## Create a General table of mesurements
	alldata <- cbind(visdata, fitdata)
	
	## Add Label columns
	alldata$Label <- rep(basename(zidir), nrow(alldata))
	## Transform Id column in !Item column
	names(alldata)[grep("Id", names(alldata))] <- "!Item"
	## Select only useful columns
	alldata$FIT_Filename <- NULL
	
	## Create _dat1.zim file
	zidatfile <- file.path(zidir, basename(zidir),
		sub("\\.zim$", "_dat1.zim", basename(zimfile)))
	if (!file.exists(dirname(zidatfile))) {
		warning("directory ", dirname(zidatfile), " does not exist")
		return(invisible(FALSE))
	}
	file.copy(from = zimfile, to = zidatfile, overwrite = FALSE)
	
	## Add table of measurements at the end
	cat("\n[Data]\n", file = zidatfile, append = TRUE)
	write.table(alldata, file = zidatfile, append = TRUE, quote = FALSE,
		sep = "\t", col.names = TRUE, row.names = FALSE, dec = ".")
	
	invisible(TRUE)
}

## Create several dat1.zim files
zimDatMakeFlowCAMAll <- function (path = ".", zimfiles = NULL)
{
	## First, switch to that directory
	if (!checkDirExists(path)) return(invisible(FALSE))
	
	if (!length(zimfiles)) { # Compute them from path
		zimfiles <- dir(path, pattern = extensionPattern("zim"),
			full.names = TRUE, all.files = FALSE, recursive = TRUE)
	}	
	## Check at least one .zim file is found or provided
	if (!length(zimfiles)) {
		warning("you must select a dir containing ZIM file(s)")
		return(invisible(FALSE))
	}
	
	## Batch process creation of all dat1.zim files
	message("Creating _dat1.zim files...")
	flush.console()
	ok <- batch(zimfiles, zimDatMakeFlowCAM, verbose = FALSE)
	if (!ok) {
		warning(sum(attr(ok, "ok")), "/", length(zimfiles),
			" files were completed (see .last.batch)")
		invisible(FALSE)
	} else {
		message("-- Done! --")
		invisible(TRUE)
	}
}

## FlowCAM special treatment because the plugin doesn't export dat1.zim!
## Read list file
## TODO: avoid duplicated code between versions
.lstRead <- function (lstfile, skip = 2)
{
	## Determine the version of the FlowCAM
	ncol <- length(read.table(lstfile, header = FALSE, sep = ":", dec = ".",
		skip = skip, nrows = 1))
	if (ncol <= 44) {
		## FlowCAM II with 44 columns
		## Read the table
		tab <- read.table(lstfile, header = FALSE, sep = ":", dec = '.',
			col.names = c("Id", "FIT_Cal_Const", "FIT_Raw_Area",
				"FIT_Raw_Feret_Max", "FIT_Raw_Feret_Min", "FIT_Raw_Feret_Mean",
				"FIT_Raw_Perim", "FIT_Raw_Convex_Perim", "FIT_Area_ABD",
				"FIT_Diameter_ABD", "FIT_Length", "FIT_Width",
				"FIT_Diameter_ESD", "FIT_Perimeter", "FIT_Convex_Perimeter",
				"FIT_Intensity", "FIT_Sigma_Intensity", "FIT_Compactness",
				"FIT_Elongation", "FIT_Sum_Intensity", "FIT_Roughness",
				"FIT_Feret_Max_Angle", "FIT_Avg_Red", "FIT_Avg_Green",
				"FIT_Avg_Blue", "FIT_PPC", "FIT_Ch1_Peak", "FIT_Ch1_TOF",
				"FIT_Ch2_Peak", "FIT_Ch2_TOF", "FIT_Ch3_Peak", "FIT_Ch3_TOF",
				"FIT_Ch4_Peak", "FIT_Ch4_TOF", "FIT_Filename", "FIT_SaveX",
				"FIT_SaveY", "FIT_PixelW", "FIT_PixelH", "FIT_CaptureX",
				"FIT_CaptureY", "FIT_High_U32", "FIT_Low_U32", "FIT_Total"),
			skip = skip)
		## Add columns present in list files from FlowCAM III
		tab$FIT_Feret_Min_Angle <- NA
		tab$FIT_Edge_Gradient <- NA
		tab$FIT_Timestamp1 <- NA
		tab$FIT_Timestamp2 <- NA
		tab$FIT_Source_Image <- NA
		tab$FIT_Calibration_Image <- NA
		tab$FIT_Ch2_Ch1_Ratio <- tab$FIT_Ch2_Peak / tab$FIT_Ch1_Peak
		## New variables calc (present in dataexport.csv from the FlowCAM)
		tab$FIT_Volume_ABD <- (4/3) * pi * (tab$FIT_Diameter_ABD/2)^3
		tab$FIT_Volume_ESD <- (4/3) * pi * (tab$FIT_Diameter_ESD/2)^3
		tab$FIT_Aspect_Ratio <- tab$FIT_Width / tab$FIT_Length
		tab$FIT_Transparency <- 1 - (tab$FIT_Diameter_ABD/tab$FIT_Diameter_ESD)
		tab$FIT_Red_Green_Ratio <- tab$FIT_Avg_Red / tab$FIT_Avg_Green
		tab$FIT_Blue_Green_Ratio <- tab$FIT_Avg_Blue / tab$FIT_Avg_Green
		tab$FIT_Red_Blue_Ratio <- tab$FIT_Avg_Red / tab$FIT_Avg_Blue
	} else { # FlowCAM III with 47 columns
		## Read the table
		tab <- read.table(lstfile, header = FALSE, sep = ":", dec = '.',
			col.names = c("Id", "FIT_Cal_Const", "FIT_Raw_Area",
				"FIT_Raw_Feret_Max", "FIT_Raw_Feret_Min", "FIT_Raw_Feret_Mean",
				"FIT_Raw_Perim", "FIT_Raw_Convex_Perim", "FIT_Area_ABD",
				"FIT_Diameter_ABD", "FIT_Length", "FIT_Width",
				"FIT_Diameter_ESD", "FIT_Perimeter", "FIT_Convex_Perimeter",
				"FIT_Intensity", "FIT_Sigma_Intensity", "FIT_Compactness",
				"FIT_Elongation", "FIT_Sum_Intensity", "FIT_Roughness",
				"FIT_Feret_Max_Angle", "FIT_Feret_Min_Angle", "FIT_Avg_Red",
				"FIT_Avg_Green", "FIT_Avg_Blue", "FIT_PPC", "FIT_Ch1_Peak",
				"FIT_Ch1_TOF", "FIT_Ch2_Peak", "FIT_Ch2_TOF", "FIT_Ch3_Peak",
				"FIT_Ch3_TOF", "FIT_Ch4_Peak", "FIT_Ch4_TOF", "FIT_Filename",
				"FIT_SaveX", "FIT_SaveY", "FIT_PixelW", "FIT_PixelH",
				"FIT_CaptureX", "FIT_CaptureY", "FIT_Edge_Gradient",
				"FIT_Timestamp1", "FIT_Timestamp2", "FIT_Source_Image",
				"FIT_Calibration_Image"),
			skip = skip)
		## Add columns present in list files from FlowCAM II
		tab$FIT_High_U32 <- NA
		tab$FIT_Low_U32 <- NA
		tab$FIT_Total <- NA
		## New variables calcul (present in dataexport.csv from the FlowCAM)
		tab$FIT_Volume_ABD <- (4/3) * pi * (tab$FIT_Diameter_ABD/2)^3
		tab$FIT_Volume_ESD <- (4/3) * pi * (tab$FIT_Diameter_ESD/2)^3
		tab$FIT_Aspect_Ratio <- tab$FIT_Width / tab$FIT_Length
		tab$FIT_Transparency <- 1 - (tab$FIT_Diameter_ABD/tab$FIT_Diameter_ESD)
		tab$FIT_Red_Green_Ratio <- tab$FIT_Avg_Red / tab$FIT_Avg_Green
		tab$FIT_Blue_Green_Ratio <- tab$FIT_Avg_Blue / tab$FIT_Avg_Green
		tab$FIT_Red_Blue_Ratio <- tab$FIT_Avg_Red / tab$FIT_Avg_Blue
		tab$FIT_Ch2_Ch1_Ratio <- tab$FIT_Ch2_Peak / tab$FIT_Ch1_Peak
	}
	tab
}

## Read context file
## TODO: avoid duplicated code between versions
.ctxRead <- function(ctxfile, fill = FALSE, largest = FALSE, vignettes = TRUE,
scalebar = TRUE, enhance = FALSE, outline = FALSE, masks = FALSE,
verbose = TRUE)
{
	## Check arguments
	if (length(ctxfile) != 1 || !is.character(ctxfile)) {
		warning("you must select one FlowCAM context (.ctx) file")
		return(NULL)
	}
	if (!checkFileExists(ctxfile, extension = "ctx", force.file = TRUE))
		return(NULL)

	## Extract information from context file
	message("reading data from FlowCAM context file '", basename(ctxfile), "'")
	## Scan the ctx file
	ctxdata <- scan(ctxfile, character(), sep = "\t", skip = 0,
		blank.lines.skip = FALSE, flush = TRUE, quiet = TRUE, comment.char = "")
	## Read version of Visual SpreadSheet
	ImageLine <- grep("^SoftwareVersion", ctxdata)
	SoftwareVersion <- as.character(sub("[ ]*$", "",
		sub("^SoftwareVersion[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
	Version <- sub("...$", "", SoftwareVersion)
	## Read right parameters
	if (Version == "1.5") {
		## Read recalibration duration
		ImageLine <- grep("^SaveIntervalMinutes", ctxdata)
		interval <- as.numeric(sub("[ ]*$", "",
			sub("^SaveIntervalMinutes[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		## Read pixel size
		ImageLine <- grep("^CalibrationConstant", ctxdata)
		pixelsize <- as.numeric(sub("[ ]*$", "",
			sub("^CalibrationConstant[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		## Read minimal size
		ImageLine <- grep("^MinESD", ctxdata)
		minsize <- as.numeric(sub("[ ]*$", "",
			sub("^MinESD[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		## Read maximal size
		ImageLine <- grep("^MaxESD", ctxdata)
		maxsize <- as.numeric(sub("[ ]*$", "",
			sub("^MaxESD[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		## Read the kind of segmentation used
		ImageLine <- grep("^CaptureDarkOrLightPixels", ctxdata)
		DarkOrLight <- as.numeric(sub("[ ]*$", "",
			sub("^CaptureDarkOrLightPixels[ ]*[=][ ]*", "",
			ctxdata[ImageLine[1]])))
		if (DarkOrLight == 0) {
			use <- "dark"
		} else if (DarkOrLight == 1) {
			use <- "light"	
		} else use <- "both"
		## Read segmentation threshold
		ImageLine <- grep("^Threshold", ctxdata)
		thresholddark <- as.numeric(sub("[ ]*$", "",
			sub("^Threshold[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		thresholdlight <- as.numeric(sub("[ ]*$", "",
			sub("^Threshold[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
    
		## Path of the export of data
		select <- file.path(basename(dirname(ctxfile)), "data_export.csv")
		## Sample name
		Sample_Name <- basename(dirname(ctxfile))
		## Read Fluo information
		ImageLine <- grep("^Ch1Gain", ctxdata)
		Gain_Fluo_Ch1 <- as.numeric(sub("[ ]*$", "",
			sub("^Ch1Gain[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		ImageLine <- grep("^Ch1Threshold", ctxdata)
		Threshold_Fluo_Ch1 <- as.numeric(sub("[ ]*$", "",
			sub("^Ch1Threshold[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		ImageLine <- grep("^Ch2Gain", ctxdata)
		Gain_Fluo_Ch2 <- as.numeric(sub("[ ]*$", "",
			sub("^Ch2Gain[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		ImageLine <- grep("^Ch2Threshold", ctxdata)
		Threshold_Fluo_Ch2 <- as.numeric(sub("[ ]*$", "",
			sub("^Ch2Threshold[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		## Read information about FlowCell
		ImageLine <- grep("^FlowCellDepth", ctxdata)
		FlowCell <- as.numeric(sub("[ ]*$", "",
			sub("^FlowCellDepth[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		## Distance to nearest
		ImageLine <- grep("^DistanceToNeighbor", ctxdata)
		Dist_To_Nearest <- as.numeric(sub("[ ]*$", "",
			sub("^DistanceToNeighbor[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		## Calculation of volume analyzed
		## Number of raw images analyzed
		ImageLine <- grep("^RawImageTotal", ctxdata)
		Raw <- as.numeric(sub("[ ]*$", "",
			sub("^RawImageTotal[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		## Area analysed (Length * Width) in pixels
		ImageLine <- grep("^AcceptableLeft", ctxdata)
		Left <- as.numeric(sub("[ ]*$", "",
			sub("^AcceptableLeft[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		ImageLine <- grep("^AcceptableRight", ctxdata)
		Right <- as.numeric(sub("[ ]*$", "",
			sub("^AcceptableRight[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		ImageLine <- grep("^AcceptableTop", ctxdata)
		Top <- as.numeric(sub("[ ]*$", "",
			sub("^AcceptableTop[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		ImageLine <- grep("^AcceptableBottom", ctxdata)
		Bottom <- as.numeric(sub("[ ]*$", "",
			sub("^AcceptableBottom[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		## Calculation of area of one image in
		## micron <= (R-L * PixelSize) * (B-T * PixelSize)
		Area <- ((Right - Left) * pixelsize) * ((Bottom - Top) * pixelsize)
		## Total volume analysed (cm^3 = ml)
		VolumeDigitized <- (Area/(10^8)) * (FlowCell/10000) * Raw
		## New fields in the ctx FlowCAM III
		Threshold_Scatter <- NA
		VolumeDigitized_VIS <- NA
		Dilution_VIS <- NA
		AutoImageRate <- NA
		FlashDuration <- NA
	} else if (Version == "2.1" || Version == "2.2") {
		## Fields not present in new version
		Gain_Fluo_Ch1 <- NA
		Gain_Fluo_Ch2 <- NA
		## Read recalibration duration
		ImageLine <- grep("^RecalibrationIntervalMinutes", ctxdata)
		interval <- as.numeric(sub("[ ]*$", "",
			sub("^RecalibrationIntervalMinutes[ ]*[=][ ]*", "",
			ctxdata[ImageLine[1]])))
		## Read pixel size
		ImageLine <- grep("^CalibrationConstant", ctxdata)
		pixelsize <- as.numeric(sub("[ ]*$", "",
			sub("^CalibrationConstant[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		## Read minimal size
		ImageLine <- grep("^MinESD", ctxdata)
		minsize <- as.numeric(sub("[ ]*$", "",
			sub("^MinESD[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		## Read maximal size
		ImageLine <- grep("^MaxESD", ctxdata)
		maxsize <- as.numeric(sub("[ ]*$", "",
			sub("^MaxESD[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		## Read the kind of segmentation used
		ImageLine <- grep("^CaptureDarkOrLightPixels", ctxdata)
		CaptureDarkOrLightPixels <- as.numeric(sub("[ ]*$", "",
			sub("^CaptureDarkOrLightPixels[ ]*[=][ ]*", "",
			ctxdata[ImageLine[1]])))
		if (CaptureDarkOrLightPixels == 0) {
			use <- "dark"
		} else if (CaptureDarkOrLightPixels == 1) {
			use <- "light"
		} else use <- "both"
		## Read segmentation threshold
		ImageLine <- grep("^ThresholdDark", ctxdata)
		thresholddark <- as.numeric(sub("[ ]*$", "",
			sub("^ThresholdDark[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		ImageLine <- grep("^ThresholdLight", ctxdata)
		thresholdlight <- as.numeric(sub("[ ]*$", "",
			sub("^ThresholdLight[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))

		## Path of the export of data
		ImageLine <- grep("^AutoExportList", ctxdata)
		AutoExportList <- as.numeric(sub("[ ]*$", "",
			sub("^AutoExportList[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		if (AutoExportList == 1) {
			select <- file.path(basename(dirname(ctxfile)),
				paste(basename(dirname(ctxfile)),"csv", sep = "."))
		} else {
			select <- file.path(basename(dirname(ctxfile)), "data_export.csv")
		}
		## Sample name
		Sample_Name <- basename(dirname(ctxfile))
		## Read Fluo information
		ImageLine <- grep("^Ch1Threshold", ctxdata)
		Threshold_Fluo_Ch1 <- as.numeric(sub("[ ]*$", "",
			sub("^Ch1Threshold[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		ImageLine <- grep("^Ch2Threshold", ctxdata)
		Threshold_Fluo_Ch2 <- as.numeric(sub("[ ]*$", "",
			sub("^Ch2Threshold[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		## Read Scatter information
		ImageLine <- grep("^ScatterThreshold", ctxdata)
		Threshold_Scatter <- as.numeric(sub("[ ]*$", "",
			sub("^ScatterThreshold[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		## Read information about FlowCell
		ImageLine <- grep("^FlowCellDepth", ctxdata)
		FlowCell <- as.numeric(sub("[ ]*$", "",
			sub("^FlowCellDepth[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		## Distance to nearest
		ImageLine <- grep("^DistanceToNeighbor", ctxdata)
		Dist_To_Nearest <- as.numeric(sub("[ ]*$", "",
			sub("^DistanceToNeighbor[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		## Calculation of volume analyzed
		## Number of raw images analyzed
		ImageLine <- grep("^RawImageTotal", ctxdata)
		Raw <- as.numeric(sub("[ ]*$", "",
			sub("^RawImageTotal[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		## Area analysed (Length * Width) in pixels
		ImageLine <- grep("^AcceptableLeft", ctxdata)
		Left <- as.numeric(sub("[ ]*$", "",
			sub("^AcceptableLeft[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		ImageLine <- grep("^AcceptableRight", ctxdata)
		Right <- as.numeric(sub("[ ]*$", "",
			sub("^AcceptableRight[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		ImageLine <- grep("^AcceptableTop", ctxdata)
		Top <- as.numeric(sub("[ ]*$", "",
			sub("^AcceptableTop[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		ImageLine <- grep("^AcceptableBottom", ctxdata)
		Bottom <- as.numeric(sub("[ ]*$", "",
			sub("^AcceptableBottom[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		## Calculation of area of one image in microns
		## <= (R-L * PixelSize) * (B-T * PixelSize)
		Area <- ((Right - Left) * pixelsize) * ((Bottom - Top) * pixelsize)
		## Total volume analysed (cm^3 = ml)
		VolumeDigitized <- (Area/(10^8)) * (FlowCell/10000) * Raw
		## Total volume analyzed calculated by Visual Spreadsheet
		ImageLine <- grep("^TotalVolumeML", ctxdata)
		VolumeDigitized_VIS <- as.numeric(sub("[ ]*$", "",
			sub("^TotalVolumeML[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		## Dilution
		ImageLine <- grep("^VolumeCalibrationFactor", ctxdata)
		Dilution_VIS <- as.numeric(sub("[ ]*$", "",
			sub("^VolumeCalibrationFactor[ ]*[=][ ]*", "",
			ctxdata[ImageLine[1]])))
		## AutoImage
		ImageLine <- grep("^AutoImageRate", ctxdata)
		AutoImageRate <- as.numeric(sub("[ ]*$", "",
		sub("^AutoImageRate[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
		## Flash Duration
		ImageLine <- grep("^FlashDuration", ctxdata)
		FlashDuration <- as.numeric(sub("[ ]*$", "",
		sub("^FlashDuration[ ]*[=][ ]*", "", ctxdata[ImageLine[1]])))
	}
	## Create the table and return it
	data.frame(select, interval, pixelsize, minsize, maxsize, use,
		thresholddark, thresholdlight, fill = fill, largest = largest,
		vignettes = vignettes, scalebar = scalebar, enhance = enhance,
		outline = outline, masks = masks, verbose = verbose, Sample_Name,
		FlowCell, Gain_Fluo_Ch1, Threshold_Fluo_Ch1, Gain_Fluo_Ch2,
		Threshold_Fluo_Ch2, Threshold_Scatter, Dist_To_Nearest, VolumeDigitized,
		VolumeDigitized_VIS, SoftwareVersion, Dilution_VIS, AutoImageRate,
		FlashDuration)
}

## Read several ctx files
.ctxReadAll <- function (path = ".", ctxfiles = NULL, fill = FALSE,
largest = FALSE, vignettes = TRUE, scalebar = TRUE, enhance = FALSE,
outline = FALSE, masks = FALSE, verbose = TRUE)
{
	## First, switch to that directory
	if (!checkDirExists(path)) return(NULL)
	
	if (!length(ctxfiles)) { # Compute them from path
		ctxfiles <- dir(path, pattern = extensionPattern("ctx"),
			full.names = TRUE, all.files = FALSE, recursive = TRUE)
	}	
	## Check at least one .ctx file is found or provided
	if (!length(ctxfiles)) {
		warning("you must select a directory containing FlowCAM data")
		return(NULL)
	}
	
	## Read first ctx file
	res <- .ctxRead(ctxfiles[1], fill = FALSE, largest = FALSE, vignettes = TRUE,
		scalebar = TRUE, enhance = FALSE, outline = FALSE, masks = FALSE,
		verbose = TRUE)
	## Make a loop to read each one
	nfiles <- length(ctxfiles)
	if (nfiles > 1) for (i in 2:nfiles)
		res <- rbind(res, .ctxRead(ctxfiles[i], fill = fill,
			largest = largest, vignettes = vignettes, scalebar = scalebar,
			enhance = enhance, outline = outline, masks = masks,
			verbose = verbose))
	res
}
