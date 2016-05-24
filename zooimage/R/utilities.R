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
## along with ZooImage. If not, see <http://www.gnu.org/licenses/>.

## Get information about a sample, given its name
sampleInfo <- function (filename,  type = c("sample", "fraction", "image",
"scs", "date", "id", "frac", "imgnbr"), ext = "_dat1[.]zim$")
{	
	base <- basename(as.character(filename))
	if (ext != "") base <- sub(ext, "", base)
	
	## Filename without extension is supposed to follow the convention:
	## scs.date.id+f[img] with scs.date.id forming an unique sample identifier
	## Note: not all verifications are conducted. So, it sometimes returns a
	## result even if the name does not conform to this specification!
	### TODO: check that the name follows the convention and determine what is
	##         optional, like date, for instance)
	switch(match.arg(type),
		sample = sub("\\+[a-zA-Z][0-9.]+$", "", base),
		fraction = sub("[0-9.]+$", "", base),
		image = base,
		scs = sub("[+.].+$", "", base),
		date = as.Date(sub("^.*([0-9]{4}-[0-1][0-9]-[0-3][0-9]).*$", "\\1", base)),
		id = sub("^.*\\..*\\.(.*)\\+.*$", "\\1", base),
		frac = sub("^.*\\+([a-zA-Z]).*$", "\\1",base),
		imgnbr = as.numeric(sub("^.*\\+[a-zA-Z]([0-9.]*)$", "\\1", base)),
		{
			warning("'type' must be 'sample', 'fraction', 'image', 'scs', 'date', 'id', 'frac' or 'imgnbr'")
			character(0)
		}
	)
}

## Convert underscores into spaces
underscoreToSpace <- function (string)
	gsub("_", " ", string)

## Trim leading and trailing white spaces and tabs
trimString <- function (string)
	sub("\\s+$", "", sub("^\\s+", "", string))

## All sample with at least one entry in a given object
listSamples <- function (ZIobj)
{ 	
	if (!inherits(ZIobj, c("ZIDat", "ZIDesc","ZITrain","ZITest"))) {
		warning("'ZIobj' must be a 'ZIDat', 'ZIDesc', 'ZITrain' or 'ZITest' object")
		return(character(0))
	}
	
	## List all samples represented in a given object
	if (inherits(ZIobj, "ZIDat")) {
    	res <- sort(unique(sampleInfo(as.character(ZIobj$Label),
			type = "sample", ext = "")))
	} else if (inherits(ZIobj, "ZIDesc")) {
		res <- sort(unique(as.character(ZIobj$Label)))
	} else if (inherits(ZIobj, c("ZITrain", "ZITest"))) {
    	res <- as.character(ZIobj$Id)
		res <- sub("_[0-9]*$", "", res)
		res <- sort(unique(sampleInfo(res, type = "sample", ext = "")))
	}
	res
}

## Unique identifiers (Ids) are a combination of Label and Item
makeId <- function (ZIDat)
	paste(ZIDat$Label, ZIDat$Item, sep = "_")

## Add classes into a ZIDat object, from ZITrain or ZITest objects
addClass <- function (ZIDat, ZIobj)
{
	## Is there a 'Class' variable in ZIobj?
	Cl <- ZIobj$Class
	if (!length(Cl))
		stop("No 'Class' column found in the ZIobj object")
	## Select only those items that are in ZIDat, in the correct order...
	Id <- ZIobj$Id
	if (!length(Id)) Id <- makeId(ZIobj)
	if (!length(Id)) stop("unable to get particle Ids from 'ZIobj'")
	names(Cl) <- Id
	ZIDat$Class <- Cl[makeId(ZIDat)]
	ZIDat
}

## Default list of variables to drop
## Version 3.0-1: added useless FIT variables
dropVars <- function ()
{
	res <- try(get("ZI.dropVarsDef"), silent = TRUE)
	if (inherits(res, "try-error"))
		res <- getOption("ZI.dropVarsDef",
			c("Id", "Label", "Item", "X", "Y", "XM", "YM", "BX", "BY", "Width",
			"Height", "Angle", "XStart", "YStart", "Dil", "Predicted",
			"Predicted2", "FIT_Cal_Const", "FIT_Avg_Red", "FIT_Avg_Green",
			"FIT_Avg_Blue", "FIT_PPC", "FIT_Ch1_Peak", "FIT_Ch1_TOF",
			"FIT_Ch2_Peak", "FIT_Ch2_TOF", "FIT_Ch3_Peak", "FIT_Ch3_TOF",
			"FIT_SaveX", "FIT_SaveY", "FIT_PixelW", "FIT_PixelH",
			"FIT_CaptureX", "FIT_CaptureY", "FIT_Edge_Gradient",    
			"FIT_Source_Image", "FIT_Calibration_Image", "FIT_High_U32",
			"FIT_Low_U32", "FIT_Total", "FIT_Red_Green_Ratio",
			"FIT_Blue_Green_Ratio", "FIT_Red_Blue_Ratio",   
			"FIT_Ch2_Ch1_Ratio", "FIT_Ch4_Peak", "FIT_Ch4_TOF", "FIT_Timestamp1",
			"FIT_Timestamp2", "FIT_Camera", "FIT_FringSize", "FIT_CircleFit",
			"FIT_Ch1_Area", "FIT_Ch2_Area", "FIT_Ch3_Area",         
			"FIT_TimeStamp1", "FIT_Source_Image.1"))
	as.character(res)
}


## Calculate derived variables... default function
calcVars <- function (x, drop.vars = NULL, drop.vars.def = dropVars())
{	
	## This is the calculation of derived variables
	## Note that you can make your own version of this function for more
	## calculated variables!
	
	## A small hack to correct some 0 (which can be problematic in further calcs)
	noZero <- function (x) {
		x[x == 0] <- 0.000000001
		x
	}
	## Euclidean distance between two points
	distance <- function (x, y)
		sqrt(x^2 + y^2)
	
	x$Minor <- noZero(x$Minor)
	x$Major <- noZero(x$Major) 
	x$AspectRatio <- x$Minor / x$Major 
	x$CentBoxD <- distance(x$BX + x$Width/2 - x$X, x$BY + x$Height/2 - x$Y)
	x$GrayCentBoxD <- distance(x$BX + x$Width/2 - x$XM, x$BY + x$Height/2 - x$YM)
	x$CentroidsD <- distance(x$X - x$XM, x$Y - x$YM)
	x$Range <- x$Max - x$Min
	x$MeanPos <- (x$Max - x$Mean) / x$Range
	x$SDNorm <- x$StdDev / x$Range
	x$CV <- x$StdDev / x$Mean * 100
	x$Area <- noZero(x$Area)
	#x$logArea <- log(x$Area)
	x$Perim. <- noZero(x$Perim.)
	#x$logPerim. <- log(x$Perim.)
	#x$logMajor <- log(x$Major)
	#x$logMinor <- log(x$Minor)
	#x$logECD <- log(noZero(x$ECD))
	x$Feret <- noZero(x$Feret)
	#x$logFeret <- log(x$Feret)
	x$MeanDia <- (x$Major + x$Minor) / 2
	x$MeanFDia <- (x$Feret + x$Minor) / 2
	#x$logMeanDia <- log(x$MeanDia)
	#x$logMeanFDia <- log(x$MeanFDia)
	x$Transp1 <- 1 - (x$ECD / x$MeanDia)
	x$Transp1[x$Transp1 < 0] <- 0
	x$Transp2 <- 1 - (x$ECD / x$MeanFDia)
	x$Transp2[x$Transp2 < 0] <- 0
	PA <- x$Perim.^2/16 - x$Area
	x$Elongation <- ifelse(PA <= 0, 1, x$Area / (x$Perim./4 - PA^.5)^2)
	x$Compactness <-  x$Perim.^2/4/pi/x$Area  # env. 1/Circ.
	x$Roundness <- 4 * x$Area / (pi * sqrt(x$Major))
	
	## Eliminate variables that are not predictors... and use Id as rownames
	Id <- x$Id
	if (length(Id)) rownames(x) <- Id
	
	## Variables to drop
	dropAll <- unique(as.character(c(drop.vars, drop.vars.def)))
	for (dropVar in dropAll) x[[dropVar]] <- NULL
	
	## Return the recalculated data frame
	x
}

## Calculate equivalent circular diameter (similar to equivalent spherical
## diameter, but for 2D images)
ecd <- function (area)
	2 * sqrt(area / pi)

## Parse an ini file (.zim, .zie, etc. are .ini files!)
### TODO: manage the case where there is no '=' in the data!
parseIni <- function (data, label = "1")
{
	## Parse an ini file (tag=value => 'tag', 'value')
	## and make a list with different sections
	
	# Is str a section?
	is.section <- function (str)
		as.logical(length(grep("^\\[.+\\]$", trimString(str)) > 0))

	## Get the name of a section
	get.section.name <- function (str)
		sub("^\\[", "", sub("\\]$", "", trimString(str)))

	## Transform a vector of characters into a data frame,
	## possibly with type conversion
	vector.convert <- function (vec)
		as.data.frame(lapply(as.list(vec), type.convert))

	if (!length(data) || !inherits(data, "character"))
		return(character(0))
	
	## Trim leading and trailing white spaces
	data <- trimString(data)
	
	## Convert underscore to space
	data <- underscoreToSpace(data)
	
	## Eliminate empty lines
	data <- data[data != ""]
	data <- paste(data, " ", sep = "")
	if (!length(data)) return(character(0))
	## Substitute the first '=' sign by another separator unlikely to appear in
	## the argument
	data <- sub("=", "&&&&&", data)
	
	## Split the strings according to this separator
	data <- strsplit(data, "&&&&&")
	
	## Get a matrix
	data <- t(as.data.frame(data))
	rownames(data) <- NULL
	
	## Make sure we have a section for the first entries (otherwise, use [.])
	if (!is.section(data[1, 1]))
		data <- rbind(c("[.]", "[.]"), data)
	Names <- as.vector(trimString(data[, 1]))
	Dat <- as.vector(trimString(data[, 2]))
	
	## Determine which is a section header
	Sec <- grep("\\[.+\\]$", Names)
	SecNames <- get.section.name(Names[Sec])
	
	## Make a vector of sections
	if (length(Sec) == 1) {
		SecNames <- rep(SecNames, length(Names))
	} else {
		SecNames <- rep(SecNames, c(Sec[2:length(Sec)],
			length(Names) + 1) - Sec)
	}
	
	## Replace section headers from all vectors
	Names[Sec] <- "Label"
	Dat[Sec] <- label
	names(Dat) <- Names
	
	## Transform SecNames in a factor
	SecNames <- as.factor(SecNames)
	
	## Split Dat on sections
	DatSec <- split(Dat, SecNames)
	
	## For each section, transform the vector in a data frame and possibly
	## convert its content
	DatSec <- lapply(DatSec, vector.convert)
	
	## Eliminate "Label" if it is ""
	if (label == "")
		DatSec <- lapply(DatSec, function(x) x[-1])
	
	DatSec	
}

## Garyscale calibration in O.D. scale
## TODO: rework all this using ImageJ
calibrate <- function (ODfile)
{
	### TODO: include also a spatial calibration procedure
	## (with a black circle around the center of the image)
	## and check also other characteristics, especially the sharpness

    cal <- c(NA, NA)
	names(cal) <- c("WhitePoint", "BlackPoint")
	msg <- character(0)

	if (!file.exists(ODfile)) {
		msg <- paste("O.D. file '", ODfile, "' not found!", sep = "")
		attr(cal, "msg") <- msg
		return(cal)
	}

	## Is it a test file?
	if (.isTestFile(ODfile)) {
		## We behave like if the file was correct and return fake calibration data!
        cal <- c(1000, 50000)
		names(cal) <- c("WhitePoint", "BlackPoint")
		attr(cal, "msg") <- character(0)
		return(cal)
	}

	filedir <- dirname(ODfile)
	if (filedir != ".") {
		## Temporary change directory to the one where the file is located
		inidir <- setwd(filedir)
		on.exit(setwd(inidir))
		ODfile <- basename(ODfile)
	}
	
	## The command to use depends on the format of the image (determined on the
	## extension)
	ext <- tolower(rev(strsplit(ODfile, "\\.")[[1]])[1])
	pgmfile <- ODfile
	if (ext == "tif") {
		## First, convert into a .pgm file
		pgmfile <- paste(ODfile, "pgm", sep = ".")
####		netpbm_tifftopnm( ODfile, pgmfile )
		delfile <- TRUE
		ext <- "pgm"
	} else delfile <- FALSE
	if (ext != "pgm")
		return(paste("Unrecognized image format for '", ODfile, "'", sep = ""))
####	OD <- netpbm_pgmhist(pgmfile, delete = delfile)
	
	## Make sure we work with 16bit images
	if (max(OD$Gray) < 256) {
		msg <- c(msg, "O.D. seems to be a 8bit image (16bit required)")	
	} else {
		## Eliminate values with low number of points
		OD <- OD[OD$Count > 100, ]
		
		## Look at range: should be widespread enough, but without saturation
		rngOD <- range(OD$Gray)
		if (rngOD[2] > 65500) msg <-
			c(msg, "Images are overexposed, or whitepoint is already calibrated")
		if (rngOD[2] < 55000)
			msg <- c(msg, "Images are underexposed")
		
		## Saturation on the left-side of the histogram is not much a problem!
		if (rngOD[2] - rngOD[1] < 40000)
			msg <- c(msg, "Images lack contrast")
		## We should end up with four segments
		graylev <- OD$Gray
		gap <- (diff(graylev) > 500)
		
		## There are not *exactly* four gaps => problem with the image!
		if (sum(gap) != 4) {
			msg <- c(msg, "Impossible to calibrate O.D.: wrong image")
		} else {
			## Get the five peaks, analyze them and get modes for blank, NDx2,
			## NDx4 and NDx8
			peaks <- as.factor(cumsum(c(0, gap)) + 1)
			peaksgray <- split(graylev, peaks)
			names(peaksgray) <- c("Black", "NDx8", "NDx4", "NDx2", "White")
			
			## These are supposed to be all narrow peaks... check this
			peakspan <- sapply(peaksgray, range)
			peaksrange <- peakspan[2, ] - peakspan[1, ]
			
			## 1.2-2: width of black peak is much larger for Epson 4990
			## => be more tolerant for that peak
			if (any(peaksrange > c(20000, rep(5000, 4)))) {
				wrongpeaks <- paste(names(peaksrange)[peaksrange > 5000],
					collapse = ", ")
				msg <- c(msg, paste("Wrong O.D. image: lack of homogeneity for",
					wrongpeaks))
			}
			
			## Look for the gray levels at the top of the peaks
			peaksheight <- split(OD$Count, peaks)
			names(peaksheight) <- c("Black", "NDx8", "NDx4", "NDx2", "White")
			findmax <- function(x) which.max(lowess(x, f = 0.05, iter = 1)$y)
			peaksval <- sapply(peaksheight, findmax)
			
			## Get the number of pixels in the white peak
			nbrwhite <- peaksheight$White[peaksval["White"]]
            
			## Replace the location by the actual gray level
			for (i in 1:5)
				peaksval[i] <- peaksgray[[i]][peaksval[i]]
			## If the number of pixels for pure white is larger than the white
			## peak found, replace it by pure white (65535)
			nbrpurewhite <- OD[nrow(OD), 2] 
			if (nbrpurewhite > nbrwhite)
				peaksval["White"] <- 65535
			
			## Now, we need to calibrate the black and white points
			WhitePoint <- 65535 - peaksval["White"]
			
			## Perform a correction for the white point
			peaksval <- peaksval + WhitePoint
			
			## Transform those gray levels into O.D.
			peaksOD <- log(peaksval) * 65535 / log(65535)
			
			## Create a data frame with gray levels and corresponding OD for
			## White, NDx2, NDx4 and NDx8
			calib <- data.frame(Gray = peaksOD[5:2], OD = c(0, 0.3, 0.6, 0.9))
			
			## Fit a line on these data
			calib.lm <- lm(OD ~ Gray, data = calib)
			
			## Check that calibration line is fine (i.e., the ANOVA should
			## reject H0 at alpha = 5%)
			if (anova(calib.lm)[["Pr(>F)"]][1] > 0.01)
				msg <- c(msg, "Wrong OD calibration: not a straight line relation at alpha level = 0.01")
			
			## Check also that R squared is at least 0.98
			rsq <- summary(calib.lm)$r.squared
			if (rsq < 0.98)
				msg <- c(msg, paste("Bad OD calibration (R squared = ",
					formatC(rsq, digits = 3), ")", sep = ""))
			
			## Check linearity of the relationship by fitting a second order
			## polynome and by looking at the t-test for the x square parameter
			calib2.lm <- lm(OD ~ I(Gray^2) + Gray, data = calib)
			if (summary(calib2.lm)$coefficients["I(Gray^2)", "Pr(>|t|)"] < 0.01)
				msg <- c(msg, "Nonlinear OD calibration at alpha level = 0.01")
			
			## Calculate the value of the black point to get 0.004 OD per gray
			## level after conversion (see the manual)
			ccoef <- coef(calib.lm)
			BlackPoint <- (1.024 - ccoef[1]) / ccoef[2]
			
			## Get the calibration data
			cal[1] <- round(WhitePoint)
			cal[2] <- round(BlackPoint)						
		}
	}
	attr(cal, "msg") <- msg
	return(cal)
}
## example:
## setwd("g:/zooplankton/madagascar2macro")
## calibrate("test.tif")

## Decimal separator to use in import/export ZooImage files
getDec <- function ()
{
	Dec <- getOption("OutDec", ".")
	## It must be either "." or ","!
	if (!Dec %in% c(".", ",")) Dec <- "."
	Dec
}

## Add a comment (from a zimfile) into a zip archive
zipNoteAdd <- function (zipfile, zimfile)
{
	zipfile <- as.character(zipfile)
	if (length(zipfile) != 1) {
		warning("exactly one 'zipfile' must be provided")
		return(FALSE)
	}
	if (!file.exists(zipfile)) {
		warning("'zipfile' not found: '", basename(zipfile), "'")
		return(FALSE)
	}
	
	zimfile <- as.character(zimfile)
	if (length(zimfile) != 1) {
		warning("exactly one 'zimfile' must be provided")
		return(FALSE)
	}
	if (!file.exists(zimfile)) {
		warning("'zimfile' not found: '", basename(zimfile), "'")
		return(FALSE)
	}
	
	if (isWin()) {
		cmd <- sprintf('%s /c type "%s" | "%s" -zq "%s" ', Sys.getenv("COMSPEC"),
			zimfile, Sys.getenv("R_ZIPCMD", "zip"), zipfile)
		res <- try(system(cmd, show.output.on.console = FALSE, invisible = TRUE,
			intern = FALSE), silent = TRUE)
	} else {
		cmd <- sprintf('zip -zq "%s" < "%s" ', zipfile, zimfile)
		res <- try(system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE,
			intern = FALSE), silent = TRUE)		
	}
	if (inherits(res, "try-error")) {
		warning(as.character(res)) # Turn error into warning
		return(FALSE)
	}
	
	if (res != 0) {
		warning("error while adding .zim data to '", basename(zipfile), "'")
		FALSE
	} else TRUE
}

## Extract the comment from the zipfile
zipNoteGet <- function (zipfile, zimfile = NULL)
{
	zipfile <- as.character(zipfile)
	if (length(zipfile) != 1) {
		warning("exactly one 'zipfile' must be provided")
		return(NULL)
	}
	if (!file.exists(zipfile)) {
		warning("'zipfile' not found: '", basename(zipfile), "'")
		return(NULL)
	}
	
	if (length(zimfile)) {
		zimfile <- as.character(zimfile)
		if (length(zimfile) != 1) {
			warning("exactly one 'zimfile' must be provided")
			return(NULL)
		}
	}
	## Make sure old data do not remain in zimfile
	unlink(zimfile)
	
	## We use unzip... and assume it is located at the same place as zip!
	if (isWin()) {
		zippgm <- Sys.getenv("R_ZIPCMD", "zip")
		unzippgm <- sub("zip$", "unzip", zippgm)
		if (unzippgm == zippgm || inherits(try(system("unzip", intern = TRUE),
			silent = TRUE), "try-error")) {
			warning("'unzip' program is required, but not found")
			return(NULL)
		}
		cmd <- sprintf('"%s" -zq "%s"', unzippgm, zipfile)
		res <- try(system(cmd, invisible = TRUE, intern = TRUE), silent = TRUE)
	} else { # Linux or Mac OS X
		cmd <- sprintf('unzip -zq "%s"', zipfile)
		res <- try(system(cmd, intern = TRUE), silent = TRUE)
	}
	if (inherits(res, "try-error")) {
		warning(as.character(res))
		return(NULL)
	}
	
	if (length(res) < 2) {
		warning("no comment data found in '", basename(zipfile), "'")
		return(character(0))
	}

	## Write the output to the file if needed and return the result
	if (length(zimfile)) {
		cat(res, file = zimfile, sep = "\n")
		invisible(res)
	} else res	
}
