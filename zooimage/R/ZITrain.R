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

template <- function (object, ...)
	UseMethod("template")
	
template.default <- function (object, ...)
	attr(object, "path")

## Get the subpath of vignettes giving their classes
.getPath <- function (x, rootdir = NULL, ext = "jpg", path, classes, ...) {
## Possibly get the classification of the particles
	if (length(classes)) {
		if (inherits(classes, "function")) {
			## Run this function for getting classes
			res <- classes(x, path, ...)
		} else if (inherits(classes, "mlearning")) {
			## Use this object to predict classes
			res <- predict(classes, x, ...)
		} else if (inherits(classes, "character")) {
			## Look for one or two columns with these names
			if (length(classes) > 2) {
				warning("cannot use more than two columns for classes... using first two only")
				classes <- classes[1:2]
			}
			nms <- names(x)
			if (any(!classes %in% nms))
				stop("classes are not existing variable names")
			res <- x[[classes[1]]]
			if (length(classes) == 2) {
				isMissing <- is.na(res)
				res[isMissing] <- x[[classes[2]]][isMissing]
			}
		}
		res <- as.character(res)
		
		## Transform this into a subpath
		subpath <- unique(c("_", path))
		while (!all(path == ".")) {
			path <- unique(dirname(path))
			subpath <- c(subpath, path)
		}
		subpath <- sort(subpath[subpath != "."])
		names(subpath) <- basename(subpath)
		res <- subpath[res]
		
		## Missing data are transformed into '_'
		res[is.na(res)] <- "_"
		names(res) <- NULL
		
	} else res <- rep("_", NROW(x)) # Default to put everything in '_'
			
	if (!length(rootdir)) rootdir <- "."
	res <- file.path(rootdir, res, paste(makeId(x), ext, sep = "."))
	res
}

## Prepare 'dir\subdir' for a manual classification by expanding all vignettes
## from a given number of zidfiles to the '_' subdir, and making
## a template for subdirs
## TODO: verify that template matches classes if classes is not NULL
prepareTrain <- function (traindir, zidbfiles,
template = c("[Basic]", "[Detailed]", "[Very detailed]"), classes = NULL, ...)
{
	## First, check that dirname of traindir is valid
	if (!checkDirExists(dirname(traindir))) return(invisible(FALSE))

	if (!checkEmptyDir(traindir,
		message = 'dir "%s" is not empty. Use AddToTrain() instead!'))
		return(invisible(FALSE))

	## Then, check that all zidfiles or zidbfiles exist
	if (hasExtension(zidbfiles[1], "zidb")) dbext <- "zidb" else dbext <- "zid"
    if (!checkFileExists(zidbfiles, dbext)) return(invisible(FALSE))
    zmax <- length(zidbfiles)

	## Also look for the template
	## If the object has a path template, use it...
	path <- attr(template, "path")
	if (!length(path)) { # Look for a .zic file with classes
		template <- as.character(template)[1]
		rx <- "^[[](.+)[]]$"
		if (grepl(rx, template)) {
			## This should be a template file in the default directory
			template <- paste(sub(rx, "\\1", template), ".zic",
				sep = "")
			template <- file.path(getTemp("ZIetc"), template)
			if (!file.exists(template)) {
				warning("The file '", template, "' is not found")
				return(invisible(FALSE))
			}
		}
		## Check that this is a .zic file
		if (!zicCheck(template)) return(invisible(FALSE))
	
		## Create the other directories
		path <- scan(template, character(), sep = "\n", skip = 2,
			quiet = TRUE)
		if (!length(path)) {
			warning(sprintf("'%s' is empty or corrupted!", template))
			return(invisible(FALSE))	
		}
	}

	## Create '_' subdir
	dir_ <- file.path(traindir, "_")
	if (!forceDirCreate(dir_)) return(invisible(FALSE))

	## Create subdirectories representing classes hierarchy as in template
	message("Making directories...")
	fullpath <- file.path(traindir, path)
	for (i in 1:length(fullpath)) {
		#message(fullpath[i])
		dir.create(fullpath[i], recursive = TRUE)
	}

	## Place the vignettes...
	message("Extracting data and vignettes ...")
	flush.console()
	for (i in 1:zmax) {
		progress(i, zmax)
		if (dbext != "zidb") {
            ## Using a temporary directory to unzip all files and then copy
    		## the RData files to the train directory
    		td <- tempfile()
    		unzip(zipfile = zidbfiles[i], exdir = td)
    		datafiles <- file.path(td, list.files(td,
    			pattern = extensionPattern(".RData"), recursive = TRUE))
			if (length(datafiles)) file.copy(datafiles, traindir)
			## Get path for the vignettes and copy them there
			zidat <- zidDatRead(zidbfiles[i])
			vigpath <- .getPath(zidat, rootdir = traindir, ext = "jpg",
				path = path, classes = classes, ...)
    		names(vigpath) <- basename(vigpath)
			vignettes <- file.path(td, list.files(td,
    			pattern = extensionPattern(".jpg"), recursive = TRUE))
			if (length(vignettes)) {
				vigpath <- vigpath[basename(vignettes)]
				isMissing <- is.na(vigpath)
				vigpath[isMissing] <- file.path(dir_,
					basename(vignettes))[isMissing]
				file.copy(vignettes, vigpath)
			} else warning("no vignettes found for ", zidbfiles[i])
    		unlink(td, recursive = TRUE)
		} else {  # Use .zidb files
            ## Link .zidb database to R objects in memory
            Zidb <- zidbLink(zidbfiles[i])
            AllItems <- ls(Zidb)
            Vigns <- AllItems[-grep("_dat1", AllItems)]
            ## Extract all vignettes in their class subdirectory
            imgext <- Zidb[[".ImageType"]]
			## Get path for the vignettes and copy them there
			zidat <- zidbDatRead(zidbfiles[i])
			vigpath <- .getPath(zidat, rootdir = traindir, ext = imgext,
				path = path, classes = classes, ...)
    		names(vigpath) <- sub(paste("\\.", imgext, "$", sep = ""), "",
				basename(vigpath))
			if (length(Vigns)) {
				vigpath <- vigpath[Vigns]
				for (j in 1:length(Vigns)) {
				    vigfile <- vigpath[i]
					if (is.na(vigfile)) vigfile <- file.path(dir_,
						paste(Vigns[i], imgext, sep = "."))
					writeBin(Zidb[[Vigns[j]]], vigpath[j])
				}
			} else warning("no vignettes found for ", zidbfiles[i])
            ## Save vignettes
            ZI.sample <- Zidb$.Data
            save(ZI.sample, file = file.path(traindir, paste(sub(".zidb", "",
				basename(zidbfiles[i])), "_dat1.RData", sep = "")))
		}
	}
	progress(101) # Clear progression indicator
	message(" -- Done! --")
	invisible(TRUE)
}

## TODO: apply selection for active learning by partial validation
prepareTest <- function (testdir, zidbfiles, template, classes = NULL, ...)
{
	if (!is.null(attr(template, "path"))) template <- attr(template, "path")
	tpl <- structure(1, path = template)
	res <- prepareTrain(testdir, zidbfiles = zidbfiles,
		template = tpl, classes = classes, ...)
	## Add a .zic file there to make sure to respect training set classes
	cat("ZI3\n[path]\n", paste(template, collapse = "\n"), "\n", sep = "",
		file = file.path(testdir, "_template.zic"))
	
	invisible(res)
}

## Function to add new vignettes in a training set
addToTrain <- function (traindir, zidbfiles, classes = NULL, ...)
{
	## Check if selected zid(b) files are already classified in the training set
	Rdata <- list.files(traindir, pattern = "[.]RData$")
	RdataNew <- paste0(sub("[.]zidb?$", "", basename(zidbfiles)), "_dat1.RData")
	NewZidb <- !RdataNew %in% Rdata
	
	if (!any(NewZidb)) { # All zidbs are already in the training set
		warning("All selected ZID(B) files already in the training set")
		return(invisible(FALSE))
	} else { # Keep only new zid(b) files
		zidbfiles <- zidbfiles[NewZidb]
		warning("You have selected ", length(zidbfiles), " new ZID(B) files.\n",
			"The others files are already included in the training set")
	}
	
	## Extract vignettes to a new subdir in '_' and .RData to parent directory
	NewDir <- "_/_NewVignettes1"
	## Check if the new directory name already exists
	if (file.exists(file.path(traindir, NewDir))) {
		DirLst <- dir(file.path(traindir, "_"), pattern = "_NewVignettes")
		NewDir <- paste("_/_NewVignettes", (length(DirLst) + 1), sep = "")
	}
	
	## Check if NewDir exist
	ToPath <- file.path(traindir, NewDir)
	if (!file.exists(ToPath))
		if (!forceDirCreate(ToPath)) return(invisible(FALSE))
	
	## Extract RData in the root directory
	zmax <- length(zidbfiles)
	message("Adding data and vignettes to the training set...\n")
	for (i in 1:zmax) {
		progress(i, zmax)
		## treatment depends if it is a .zid or .zidb file
		zidbfile <- zidbfiles[i]
		if (grepl("[.]zidb$", zidbfile)) { # .zidb file
			## Link .zidb database to R objects in memory
            Zidb <- zidbLink(zidbfile)
            AllItems <- ls(Zidb)
            Vigns <- AllItems[-grep("_dat1", AllItems)]
            ## Copy all vignettes in the TopPath directory
            imgext <- Zidb[[".ImageType"]]
			## Get path for the vignettes and copy them there
			zidat <- zidbDatRead(zidbfile)
			vigpath <- .getPath(zidat, rootdir = traindir, ext = imgext,
				path = attr(zidat, "path"), classes = classes, ...)
    		vigpath[vigpath == "_"] <- ToPath
			names(vigpath) <- sub(paste("\\.", imgext, "$", sep = ""), "",
				basename(vigpath))
			if (length(Vigns)) {
				vigpath <- vigpath[Vigns]
				for (j in 1:length(Vigns)) {
				    vigfile <- vigpath[i]
					if (is.na(vigfile)) vigfile <- file.path(ToPath,
						paste(Vigns[i], imgext, sep = "."))
					writeBin(Zidb[[Vigns[j]]], vigpath[j])
				}
			} else warning("no vignettes found for ", zidbfile)
            ## Save RData file
            ZI.sample <- Zidb$.Data
            save(ZI.sample, file = file.path(traindir, paste(sub(".zidb", "",
				basename(zidbfile)), "_dat1.RData", sep = "")))

		} else { # .zid file
			## Using a temporary directory to unzip all files and then copy
			## the RData files to the train directory
			td <- tempfile()
			unzip(zipfile = zidbfile, exdir = td)
			datafiles <- file.path(td, list.files(td,
				pattern = extensionPattern(".RData"), recursive = TRUE))
			if (length(datafiles))
				file.copy(datafiles, file.path(traindir, basename(datafiles)))
			## Get path for the vignettes and copy them there
			zidat <- zidDatRead(zidbfile)
			vigpath <- .getPath(zidat, rootdir = traindir, ext = "jpg",
				path = attr(zidat, "path"), classes = classes, ...)
    		vigpath[vigpath == "_"] <- ToPath
			names(vigpath) <- basename(vigpath)
			vignettes <- file.path(td, list.files(td,
    			pattern = extensionPattern(".jpg"), recursive = TRUE))
			if (length(vignettes)) {
				vigpath <- vigpath[basename(vignettes)]
				isMissing <- is.na(vigpath)
				vigpath[isMissing] <- file.path(ToPath,
					basename(vignettes))[isMissing]
				file.copy(vignettes, vigpath)
			} else warning("no vignettes found for ", zidbfile)
			unlink(td, recursive = TRUE)	
		}
	}
	progress(101) # Clear progression indicator
	message("-- Done --\n")
	invisible(TRUE)
}

addToTest <- function (testdir, zidbfiles, classes = NULL, ...)
	invisible(addToTrain(traindir = testdir, zidbfiles = zidbfiles,
		classes = classes, ...))

## Retrieve information from a manual training set in a 'ZITrain' object	
## TODO: check dir names are unique, check no duplicated vignettes,
##       check all measurements are there, ... + exhaustive report!
getTrain <- function (traindir, creator = NULL, desc = NULL, keep_ = FALSE,
na.rm = FALSE)
{
	## 'traindir' must be the base directory of the manual classification
	if (!checkDirExists(traindir)) return(invisible(FALSE))

	## Make sure we have .RData files in this traindir (otherwise it is
	## perhaps not a training set root dir!
	Dats <- list.files(traindir, pattern = "_dat1[.]RData$", full.names = TRUE)
	if (!length(Dats)) {
		warning("'traindir' does not appear to be a ", getTemp("ZIname"),
			" training set root dir!")
		return(invisible(FALSE))
	}

	## List the .jpg or .png files (recursively) in the dir
	res <- jpgList(traindir, recursive = TRUE)
	if (!length(res)) res <- pngList(traindir, recursive = TRUE)

	## Check the result...
	if (!length(res)) {
		warning("no PNG or JPEG vignettes found in this tree")
		return(invisible(FALSE))
	}

	## Replace "\\" by "/"
	res <- gsub("[\\]", "/", res)

	## Do we eliminate the '_' directory?
	if (!is.na(keep_) && !isTRUE(as.logical(keep_)))
		res <- grep("^[^_]", res, value = TRUE)

	## 'Id' is the name of the vignettes, minus the extension
	Id <- noExtension(res)

	## 'Path' is the directory path
	Path <- dirname(res)

	## 'Class' is the last directory where the vignettes are located
	Class <- basename(Path)
	
	## For all items in _ or one of its subdirectories, replace Class by NA
	if (is.na(keep_)) Class[grepl("^[_]", res)] <- NA

	## Create a  data frame with Id and Class
	df <- data.frame(Id = Id, Class = Class)
	df$Id <- as.character(df$Id)
	nitems <- nrow(df)

	## Read in all the .RData files from the root directory and merge them
	## Get measurement infos
	ZI.sample <- NULL
	load(Dats[1])
	Dat <- ZI.sample
	Classes <- class(Dat)
	
	## Modif Kev to free memory
	Dat <- cbind(Id = makeId(Dat), Dat)
	Dat <- merge(Dat, df, by = "Id")

	if (length(Dats) > 1) {
		for (i in 2:length(Dats)) {
			load(Dats[i])
			ZI.sample <- cbind(Id = makeId(ZI.sample), ZI.sample)
			ZI.sample <- merge(ZI.sample, df, by = "Id")
			Dat <- merge(Dat, ZI.sample, all = TRUE)
			Dat$X.Item.1 <- NULL
		}
	}
	rownames(Dat) <- 1:nrow(Dat)

	## Rename Dat in df
	df <- Dat
	## Problem if there is no remaining row in the data frame
	if (nrow(df) == 0) {
		warning("No valid item found (no vignettes with valid measurement data)")
		return(invisible(FALSE))
	}

	## Check that all items have associated measurements
	if (nrow(df) < nitems)
		warning(nitems - nrow(df),
			" vignettes without measurement data are eliminated (",
			nrow(df), " items remain in the object)")

	if (any(is.na(df)))
		if (isTRUE(as.logical(na.rm))) {
  	  		message("NAs found in the table of measurements and deleted")
  	  		df <- na.omit(df)
		} else message("NAs found in the table of measurements and left there")
	
	## Add attributes
	attr(df, "traindir") <- dir
	attr(df, "path") <- unique(Path)
	if (length(creator)) attr(df, "creator") <- creator
	if (length(desc)) attr(df, "desc") <- desc
	Classes <- c("ZI3Train", "ZITrain", Classes)
	class(df) <- Classes
	
	## Be sure that variables are numeric (sometimes, wrong importation)
#	as.numeric.Vars <- function (ZIDat, numvars) {
#	    if (is.null(numvars)) # Default values
#	        numvars <- c("ECD",
#	            "FIT_Area_ABD", "FIT_Diameter_ABD", "FIT_Volume_ABD",
#				"FIT_Diameter_ESD", "FIT_Volume_ESD", "FIT_Length", "FIT_Width",
#				"FIT_Aspect_Ratio", "FIT_Transparency", "FIT_Intensity",
#				"FIT_Sigma_Intensity", "FIT_Sum_Intensity", "FIT_Compactness",
#				"FIT_Elongation", "FIT_Perimeter", "FIT_Convex_Perimeter",
#				"FIT_Roughness", "FIT_Feret_Max_Angle", "FIT_PPC", "FIT_Ch1_Peak",
#				"FIT_Ch1_TOF", "FIT_Ch2_Peak", "FIT_Ch2_TOF", "FIT_Ch3_Peak",
#				"FIT_Ch3_TOF", "FIT_Avg_Red", "FIT_Avg_Green", "FIT_Avg_Blue",
#				"FIT_Red_Green_Ratio", "FIT_Blue_Green", "FIT_Red_Blue_Ratio",
#				"FIT_CaptureX", "FIT_CaptureY", "FIT_SaveX", "FIT_SaveY",
#				"FIT_PixelW", "FIT_PixelH", "FIT_Cal_Const",
#	            "Area", "Mean", "StdDev", "Mode", "Min", "Max", "X", "Y", "XM",
#	            "YM", "Perim.", "BX", "BY", "Width", "Height", "Major", "Minor",
#				"Angle", "Circ.", "Feret", "IntDen", "Median", "Skew", "Kurt",
#				"XStart", "YStart", "Dil")
#
#	    ## Make sure numvars are numeric
#		Names <- names(ZIDat)
#	    for (numvar in numvars) {
#	        if (numvar %in% Names && !is.numeric(ZIDat[, numvar]))
#	            ZIDat[, numvar] <- as.numeric(ZIDat[, numvar])
#	    }
#	    ZIDat
#	}
#	as.numeric.Vars(df, numvars = numvars)

	df
}

getTest <- function (testdir, creator = NULL, desc = NULL, keep_ = NA,
na.rm = FALSE)
{
	## Same as getTrain() but returns a ZITest object... and read _template.zic
	## to make sure that path and classes do match!
	zicfile <- file.path(testdir, "_template.zic")
	if (!file.exists(zicfile) || !zicCheck(zicfile))
		stop("testdir does not seem to contain a valid test set (may be use getTrain()?)")
	
	res <- getTrain(traindir = testdir, creator = creator, desc = desc,
		keep_ = keep_, na.rm = na.rm)
	class(res) <- c("ZI3Test", "ZITest", class(res)[-(1:2)])
	
	## Read the _template.zic file and change res$Class factors and path accordingly
	path <- scan(zicfile, character(), sep = "\n", skip = 2, quiet = TRUE)
	if (!length(path))
		stop(sprintf("'%s' is empty or corrupted!", zicfile))
	attr(res, "path") <- path
	
	## Now, make sure to recode res$Class factor in the correct order!
	lev <- sort(basename(path))
	res$Class <- factor(as.character(res$Class), levels = lev, exclude = NULL)

	invisible(res)
}

.recodeLevels <- function (object, depth = 1)
{
	if (!inherits(object, "ZITrain"))
		stop("'ZITrain' must be a 'ZITrain' object")
	
	depth <- as.integer(depth)[1]
	
	## Get the "path" attribute
	path <- attr(object, "path")
	
	## Split strings on "/"
	path <- strsplit(path, "/", fixed = TRUE)
	
	## Functions to get last item, or an item at a given level
	level <- function (x, depth = 1)
		ifelse(length(x) >= depth, x[depth], x[length(x)])
	
	## Return a list with new levels
	sapply(path, level, depth = depth)
}

recode <- function (object, ...)
	UseMethod("recode")

recode.ZITrain <- function (object, new.levels, depth, ...)
{	
	if (!missing(depth)) {
		if (!missing(new.levels))
			warning("depth is provided, so, new.levels is ignored and recomputed")
		new.levels <- .recodeLevels(object, depth)
	}
	
	## Check that new.levels is of the same length as levels(object$Class)
	## [and object$Predicted or Predicted2, possibly]
	levels <- levels(object$Class)
	new.levels <- as.character(new.levels)
	if (length(new.levels) != length(levels))
		stop("length of new.levels must match levels in object$Class")
	
	relevel <- function (x, levels, new.levels) {
		x <- as.character(x)
		for (i in 1:length(levels))
			if (new.levels[i] != levels[i])
				x[x == levels[i]] <- new.levels[i]
		as.factor(x)
	}
	
	object$Class <- relevel(object$Class, levels, new.levels)
	if (!is.null(object$Predicted))
		object$Predicted <- relevel(object$Predicted, levels, new.levels)
	if (!is.null(object$Predicted2))
		object$Predicted2 <- relevel(object$Predicted2, levels, new.levels)
	
	## If a new path is given for these new classes, change it
	path <- attr(new.levels, "path")
	### TODO: check its validity here
	if (!is.null(path)) attr(object, "path") <- path
	object
}

recode.ZITest <- recode.ZITrain
